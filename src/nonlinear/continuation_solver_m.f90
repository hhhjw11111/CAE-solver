module continuation_solver_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE,F_PHI,global_dof
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, &
    NonlinearOptionsType, ContinuationOptionsType, DiagnosticsType
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INVALID_ARGUMENT, STATUS_CONTINUATION_FAILED
  use sparse_pattern_m, only : SparseAssemblyContextType
  use dirichlet_m, only : DirichletBCType
  use newton_solver_m, only : NewtonReportType, solve_single_load_newton
  use balance_audit_m, only : audit_global_balances
  implicit none
  private

  integer, parameter, public :: CONTINUATION_PHASE_START=0
  integer, parameter, public :: CONTINUATION_PHASE_VOLTAGE=1
  integer, parameter, public :: CONTINUATION_PHASE_JOULE=2

  type, public :: ContinuationReportType
    logical :: converged=.false.
    logical :: arclength_fold_detected=.false.
    logical :: arclength_fold_located=.false.,arclength_target_crossed=.false.
    real(dp) :: arclength_fold_bracket_width=huge(1.0_dp)
    real(dp) :: arclength_fold_left_voltage_factor=0.0_dp,arclength_fold_right_voltage_factor=0.0_dp
    real(dp) :: arclength_fold_left_tangent=0.0_dp,arclength_fold_right_tangent=0.0_dp
    real(dp) :: arclength_fold_voltage_factor=0.0_dp
    real(dp) :: arclength_fold_joule_factor=0.0_dp
    integer :: accepted_points=0
    integer :: rejected_attempts=0
    integer :: total_attempts=0
    integer :: total_newton_iterations=0
    real(dp) :: minimum_accepted_step=1.0_dp
    real(dp) :: maximum_checkpoint_change_on_rejection=0.0_dp
    real(dp) :: maximum_balance_defect=0.0_dp
    real(dp) :: maximum_boundary_error=0.0_dp
    type(LoadStateType) :: final_load
    type(LoadStateType) :: transition_load
    type(StateType) :: transition_state
    integer, allocatable :: phase(:)
    integer, allocatable :: newton_iterations(:)
    real(dp), allocatable :: path_parameter(:)
    real(dp), allocatable :: voltage_factor(:)
    real(dp), allocatable :: joule_factor(:)
    real(dp), allocatable :: residual_rms(:)
    real(dp), allocatable :: update_rms(:)
    real(dp), allocatable :: balance(:,:)
    real(dp), allocatable :: heat_free_sum(:),heat_absolute_sum(:)
    real(dp), allocatable :: boundary_error(:)
    real(dp), allocatable :: minimum_temperature(:)
    real(dp), allocatable :: maximum_temperature(:)
    real(dp), allocatable :: minimum_n(:)
    real(dp), allocatable :: maximum_n(:)
    real(dp), allocatable :: minimum_p(:)
    real(dp), allocatable :: maximum_p(:)
  end type ContinuationReportType

  public :: solve_load_continuation

  abstract interface
    subroutine accepted_checkpoint(state,load,point,residual,update,status)
      import :: StateType,LoadStateType,SolverStatusType,dp
      type(StateType), intent(in) :: state
      type(LoadStateType), intent(in) :: load
      integer, intent(in) :: point
      real(dp), intent(in) :: residual,update
      type(SolverStatusType), intent(out) :: status
    end subroutine accepted_checkpoint
  end interface

contains

  subroutine solve_load_continuation(mesh,materials,context,boundary_conditions,field_scales, &
                                     nonlinear_options,continuation_options,target_load, &
                                     state,current_load,report,status,checkpoint)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    type(NonlinearOptionsType), intent(in) :: nonlinear_options
    type(ContinuationOptionsType), intent(in) :: continuation_options
    type(LoadStateType), intent(in) :: target_load
    type(StateType), intent(inout) :: state
    type(LoadStateType), intent(inout) :: current_load
    type(ContinuationReportType), intent(out) :: report
    type(SolverStatusType), intent(out) :: status
    procedure(accepted_checkpoint), optional :: checkpoint
    type(StateType) :: trial_state,checkpoint_state
    type(StateType) :: previous_state
    type(LoadStateType) :: trial_load
    type(NonlinearOptionsType) :: attempt_options
    type(NewtonReportType) :: newton_report
    type(SolverStatusType) :: local_status
    real(dp) :: step,accepted_step,next_value,remaining,checkpoint_change
    real(dp) :: previous_parameter,previous_step,next_checkpoint
    logical :: predictor_available
    integer :: phase

    call clear_status(status)
    report=ContinuationReportType()
    call validate_continuation_inputs(field_scales,nonlinear_options, &
      continuation_options,current_load,target_load,status)
    if (.not. status%is_ok()) return
    call allocate_report(report,continuation_options%maximum_attempts+3)
    allocate(trial_state%x(size(state%x)),checkpoint_state%x(size(state%x)))
    if(continuation_options%secant_predictor) allocate(previous_state%x(size(state%x)))
    predictor_available=.false.;previous_parameter=0.0_dp
    attempt_options=nonlinear_options
    attempt_options%max_newton=min(nonlinear_options%max_newton, &
      continuation_options%maximum_newton_per_attempt)

    call solve_single_load_newton(mesh,materials,current_load,context,boundary_conditions, &
      field_scales,attempt_options,state,newton_report,local_status)
    report%total_attempts=1
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_CONTINUATION_FAILED, &
        'Continuation start state is not converged: '//trim(local_status%message))
      return
    end if
    call record_point(CONTINUATION_PHASE_START,0.0_dp,mesh,materials,context, &
      boundary_conditions,state,current_load,newton_report,report,local_status, &
      continuation_options%heat_balance_reference_power)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if

    if(present(checkpoint)) then
      call checkpoint(state,current_load,report%accepted_points,newton_report%final_residual_rms, &
        newton_report%final_update_rms,status)
      if(.not.status%is_ok()) return
    end if

    if (current_load%voltage_factor >= target_load%voltage_factor .and. &
        current_load%joule_factor == 0.0_dp) then
      call save_transition_checkpoint(state,current_load,report)
    end if

    step=continuation_options%initial_step
    phase=CONTINUATION_PHASE_VOLTAGE
    do while(current_load%voltage_factor < target_load%voltage_factor)
      call check_attempt_limit(report,continuation_options,status)
      if (.not. status%is_ok()) return
      remaining=target_load%voltage_factor-current_load%voltage_factor
      accepted_step=min(step,remaining)
      if(continuation_options%voltage_checkpoint_spacing>0.0_dp) then
        next_checkpoint=continuation_options%voltage_checkpoint_spacing* &
          real(1+floor((current_load%voltage_factor+8.0_dp*epsilon(1.0_dp))/ &
          continuation_options%voltage_checkpoint_spacing),dp)
        accepted_step=min(accepted_step,next_checkpoint-current_load%voltage_factor)
      end if
      next_value=current_load%voltage_factor+accepted_step
      trial_load=current_load
      trial_load%voltage_factor=next_value
      checkpoint_state%x=state%x
      trial_state%x=checkpoint_state%x
      if(continuation_options%secant_predictor .and.predictor_available) then
        previous_step=current_load%voltage_factor-previous_parameter
        trial_state%x=state%x+(accepted_step/previous_step)*(state%x-previous_state%x)
      else if(continuation_options%secant_predictor) then
        call lift_voltage_guess(mesh,boundary_conditions,accepted_step,trial_state)
      end if
      call solve_single_load_newton(mesh,materials,trial_load,context,boundary_conditions, &
        field_scales,attempt_options,trial_state,newton_report,local_status)
      report%total_attempts=report%total_attempts+1
      if (local_status%is_ok()) then
        if(continuation_options%secant_predictor) then
          previous_state%x=state%x;previous_parameter=current_load%voltage_factor;predictor_available=.true.
        end if
        state%x=trial_state%x
        current_load=trial_load
        report%minimum_accepted_step=min(report%minimum_accepted_step,accepted_step)
        call record_point(phase,current_load%voltage_factor,mesh,materials,context, &
          boundary_conditions,state,current_load,newton_report,report,local_status, &
          continuation_options%heat_balance_reference_power)
        if (.not. local_status%is_ok()) then
          call copy_failure(local_status,status); return
        end if
        if(present(checkpoint)) then
          call checkpoint(state,current_load,report%accepted_points,newton_report%final_residual_rms, &
            newton_report%final_update_rms,status)
          if(.not.status%is_ok()) return
        end if
        step=adapt_step(step,newton_report%iterations,continuation_options)
      else
        checkpoint_change=maxval(abs(state%x-checkpoint_state%x))
        report%maximum_checkpoint_change_on_rejection=max( &
          report%maximum_checkpoint_change_on_rejection,checkpoint_change)
        report%rejected_attempts=report%rejected_attempts+1
        step=step*continuation_options%cutback_factor
        if (step < continuation_options%minimum_step) then
          call set_status(status,STATUS_CONTINUATION_FAILED, &
            'Voltage continuation reached the minimum step after Newton failure; '// &
            'last Newton status code='//integer_to_string(local_status%code)//': '// &
            trim(local_status%message))
          return
        end if
      end if
    end do
    if (current_load%joule_factor == 0.0_dp) then
      call save_transition_checkpoint(state,current_load,report)
    end if

    step=continuation_options%initial_step
    phase=CONTINUATION_PHASE_JOULE
    predictor_available=.false.
    do while(current_load%joule_factor < target_load%joule_factor)
      call check_attempt_limit(report,continuation_options,status)
      if (.not. status%is_ok()) return
      remaining=target_load%joule_factor-current_load%joule_factor
      accepted_step=min(step,remaining)
      next_value=current_load%joule_factor+accepted_step
      trial_load=current_load
      trial_load%joule_factor=next_value
      checkpoint_state%x=state%x
      trial_state%x=checkpoint_state%x
      if(continuation_options%secant_predictor .and.predictor_available) then
        previous_step=current_load%joule_factor-previous_parameter
        trial_state%x=state%x+(accepted_step/previous_step)*(state%x-previous_state%x)
      end if
      call solve_single_load_newton(mesh,materials,trial_load,context,boundary_conditions, &
        field_scales,attempt_options,trial_state,newton_report,local_status)
      report%total_attempts=report%total_attempts+1
      if (local_status%is_ok()) then
        if(continuation_options%secant_predictor) then
          previous_state%x=state%x;previous_parameter=current_load%joule_factor;predictor_available=.true.
        end if
        state%x=trial_state%x
        current_load=trial_load
        report%minimum_accepted_step=min(report%minimum_accepted_step,accepted_step)
        call record_point(phase,current_load%joule_factor,mesh,materials,context, &
          boundary_conditions,state,current_load,newton_report,report,local_status, &
          continuation_options%heat_balance_reference_power)
        if (.not. local_status%is_ok()) then
          call copy_failure(local_status,status); return
        end if
        if(present(checkpoint)) then
          call checkpoint(state,current_load,report%accepted_points,newton_report%final_residual_rms, &
            newton_report%final_update_rms,status)
          if(.not.status%is_ok()) return
        end if
        step=adapt_step(step,newton_report%iterations,continuation_options)
      else
        checkpoint_change=maxval(abs(state%x-checkpoint_state%x))
        report%maximum_checkpoint_change_on_rejection=max( &
          report%maximum_checkpoint_change_on_rejection,checkpoint_change)
        report%rejected_attempts=report%rejected_attempts+1
        step=step*continuation_options%cutback_factor
        if (step < continuation_options%minimum_step) then
          call set_status(status,STATUS_CONTINUATION_FAILED, &
            'Joule continuation reached the minimum step after Newton failure; '// &
            'last Newton status code='//integer_to_string(local_status%code)//': '// &
            trim(local_status%message))
          return
        end if
      end if
    end do

    report%converged=.true.
    report%final_load=current_load
  end subroutine solve_load_continuation

  subroutine validate_continuation_inputs(field_scales,nonlinear_options,options, &
                                          current_load,target_load,status)
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    type(NonlinearOptionsType), intent(in) :: nonlinear_options
    type(ContinuationOptionsType), intent(in) :: options
    type(LoadStateType), intent(in) :: current_load,target_load
    type(SolverStatusType), intent(out) :: status

    call clear_status(status)
    if (any(field_scales <= 0.0_dp) .or. .not. all(ieee_is_finite(field_scales))) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Continuation field scales are invalid.')
    else if (current_load%voltage_factor < 0.0_dp .or. current_load%joule_factor < 0.0_dp .or. &
             target_load%voltage_factor > 1.0_dp .or. target_load%joule_factor > 1.0_dp .or. &
             current_load%voltage_factor > target_load%voltage_factor .or. &
             current_load%joule_factor > target_load%joule_factor .or. &
             (current_load%voltage_factor < target_load%voltage_factor .and. &
              ((.not. options%full_joule_during_voltage .and. &
                current_load%joule_factor > 0.0_dp) .or. &
               (options%full_joule_during_voltage .and. &
                current_load%joule_factor /= target_load%joule_factor)))) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Continuation load path is invalid.')
    else if (.not.ieee_is_finite(options%voltage_checkpoint_spacing) .or. &
             options%voltage_checkpoint_spacing<0.0_dp .or.options%voltage_checkpoint_spacing>1.0_dp .or. &
             (options%voltage_checkpoint_spacing>0.0_dp .and.options%voltage_checkpoint_spacing<1.0e-8_dp)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Continuation checkpoint spacing must be zero or in [1e-8,1].')
    else if (options%initial_step <= 0.0_dp .or. options%minimum_step <= 0.0_dp .or. &
             options%maximum_step < options%minimum_step .or. &
             options%initial_step > options%maximum_step .or. &
             options%cutback_factor <= 0.0_dp .or. options%cutback_factor >= 1.0_dp .or. &
             options%growth_factor <= 1.0_dp .or. options%maximum_attempts < 1 .or. &
             options%maximum_newton_per_attempt < 1 .or. &
             options%fast_iteration_threshold < 0 .or. &
             options%slow_iteration_threshold <= options%fast_iteration_threshold .or. &
             nonlinear_options%max_newton < 1) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Continuation options are invalid.')
    end if
  end subroutine validate_continuation_inputs

  subroutine allocate_report(report,capacity)
    type(ContinuationReportType), intent(inout) :: report
    integer, intent(in) :: capacity

    allocate(report%heat_free_sum(capacity),report%heat_absolute_sum(capacity))
    report%heat_free_sum=0.0_dp;report%heat_absolute_sum=0.0_dp
    allocate(report%phase(capacity),report%newton_iterations(capacity), &
      report%path_parameter(capacity),report%voltage_factor(capacity), &
      report%joule_factor(capacity),report%residual_rms(capacity), &
      report%update_rms(capacity),report%balance(4,capacity), &
      report%boundary_error(capacity),report%minimum_temperature(capacity), &
      report%maximum_temperature(capacity),report%minimum_n(capacity), &
      report%maximum_n(capacity),report%minimum_p(capacity),report%maximum_p(capacity))
    report%phase=-1
    report%newton_iterations=-1
    report%path_parameter=-1.0_dp
    report%voltage_factor=-1.0_dp
    report%joule_factor=-1.0_dp
    report%residual_rms=-1.0_dp
    report%update_rms=-1.0_dp
    report%balance=-1.0_dp
    report%boundary_error=-1.0_dp
    report%minimum_temperature=-1.0_dp
    report%maximum_temperature=-1.0_dp
    report%minimum_n=-1.0_dp
    report%maximum_n=-1.0_dp
    report%minimum_p=-1.0_dp
    report%maximum_p=-1.0_dp
  end subroutine allocate_report

  subroutine record_point(phase,path_parameter,mesh,materials,context,boundary_conditions, &
                          state,load,newton_report,report,status,heat_reference_power)
    integer, intent(in) :: phase
    real(dp), intent(in) :: path_parameter
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(NewtonReportType), intent(in) :: newton_report
    type(ContinuationReportType), intent(inout) :: report
    type(SolverStatusType), intent(out) :: status
    type(DiagnosticsType) :: diagnostics
    real(dp) :: balance(4),boundary_error
    real(dp), intent(in) :: heat_reference_power
    real(dp) :: heat_free_sum,heat_absolute_sum
    integer :: point

    call audit_global_balances(mesh,materials,state,load,context,boundary_conditions, &
      balance,boundary_error,diagnostics,status,heat_reference_power,heat_free_sum,heat_absolute_sum)
    if (.not. status%is_ok()) return
    point=report%accepted_points+1
    report%accepted_points=point
    report%phase(point)=phase
    report%path_parameter(point)=path_parameter
    report%voltage_factor(point)=load%voltage_factor
    report%joule_factor(point)=load%joule_factor
    report%newton_iterations(point)=newton_report%iterations
    report%residual_rms(point)=newton_report%final_residual_rms
    report%update_rms(point)=newton_report%final_update_rms
    report%balance(:,point)=balance
    report%heat_free_sum(point)=heat_free_sum
    report%heat_absolute_sum(point)=heat_absolute_sum
    report%boundary_error(point)=boundary_error
    report%minimum_temperature(point)=diagnostics%min_temperature
    report%maximum_temperature(point)=diagnostics%max_temperature
    report%minimum_n(point)=diagnostics%min_n
    report%maximum_n(point)=diagnostics%max_n
    report%minimum_p(point)=diagnostics%min_p
    report%maximum_p(point)=diagnostics%max_p
    report%total_newton_iterations=report%total_newton_iterations+newton_report%iterations
    report%maximum_balance_defect=max(report%maximum_balance_defect,maxval(balance))
    report%maximum_boundary_error=max(report%maximum_boundary_error,boundary_error)
  end subroutine record_point

  subroutine save_transition_checkpoint(state,load,report)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(ContinuationReportType), intent(inout) :: report

    if (allocated(report%transition_state%x)) deallocate(report%transition_state%x)
    allocate(report%transition_state%x(size(state%x)))
    report%transition_state%x=state%x
    report%transition_load=load
  end subroutine save_transition_checkpoint

  subroutine lift_voltage_guess(mesh,bc,step,state)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), intent(in) :: bc(:)
    real(dp), intent(in) :: step
    type(StateType), intent(inout) :: state
    real(dp) :: xmin,xmax,tolerance,left,right,weight
    integer :: b,node,nleft,nright
    ! Initial guess only: avoid injecting the entire voltage increment into the
    ! first contact element. Apply only to the recognized two planar x contacts.
    ! No equations, boundary values or acceptance thresholds are changed.
    xmin=minval(mesh%coordinates(1,:));xmax=maxval(mesh%coordinates(1,:))
    if(xmax<=xmin) return
    tolerance=1.0e-10_dp*(xmax-xmin);nleft=0;nright=0;left=0.0_dp;right=0.0_dp
    do b=1,size(bc)
      if(bc(b)%field_id/=F_PHI) cycle
      if(.not.allocated(bc(b)%node_ids)) return
      if(all(abs(mesh%coordinates(1,bc(b)%node_ids)-xmin)<tolerance)) then
        left=bc(b)%value_end-bc(b)%value_start;nleft=nleft+1
      else if(all(abs(mesh%coordinates(1,bc(b)%node_ids)-xmax)<tolerance)) then
        right=bc(b)%value_end-bc(b)%value_start;nright=nright+1
      else
        return
      end if
    end do
    if(nleft/=1 .or.nright/=1) return
    do node=1,mesh%number_of_nodes()
      weight=(mesh%coordinates(1,node)-xmin)/(xmax-xmin)
      state%x(global_dof(node,F_PHI))=state%x(global_dof(node,F_PHI))+step*((1.0_dp-weight)*left+weight*right)
    end do
  end subroutine lift_voltage_guess

  subroutine check_attempt_limit(report,options,status)
    type(ContinuationReportType), intent(in) :: report
    type(ContinuationOptionsType), intent(in) :: options
    type(SolverStatusType), intent(out) :: status

    call clear_status(status)
    if (report%total_attempts >= options%maximum_attempts) then
      call set_status(status,STATUS_CONTINUATION_FAILED, &
        'Continuation reached the frozen total-attempt limit.')
    end if
  end subroutine check_attempt_limit

  real(dp) function adapt_step(current_step,newton_iterations,options) result(next_step)
    real(dp), intent(in) :: current_step
    integer, intent(in) :: newton_iterations
    type(ContinuationOptionsType), intent(in) :: options

    next_step=current_step
    if (newton_iterations <= options%fast_iteration_threshold) then
      next_step=min(options%maximum_step,current_step*options%growth_factor)
    else if (newton_iterations >= options%slow_iteration_threshold) then
      next_step=max(options%minimum_step,current_step*options%cutback_factor)
    end if
  end function adapt_step

  subroutine copy_failure(source,destination)
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,source%code,trim(source%message))
  end subroutine copy_failure

  pure function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

end module continuation_solver_m

module newton_solver_m
  use, intrinsic :: iso_fortran_env, only : output_unit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE, F_ETA_N, F_ETA_P
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, &
    NonlinearOptionsType, DiagnosticsType
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INVALID_ARGUMENT, STATUS_NONLINEAR_SOLVER_FAILED
  use sparse_pattern_m, only : SparseAssemblyContextType, csr_matvec
  use dirichlet_m, only : DirichletBCType, project_dirichlet_state, build_fixed_dof_mask
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, &
    assemble_global_system
  use scaling_m, only : build_variable_scaling, equilibrate_csr_system, unscale_solution
  use pardiso_solver_m, only : PardisoDiagnosticsType, solve_csr_pardiso
  implicit none
  private

  type, public :: NewtonReportType
    logical :: converged = .false.
    integer :: iterations = 0
    integer :: history_count = 0
    integer :: total_linear_solves = 0
    integer :: total_backtracks = 0
    integer :: maximum_perturbed_pivots = 0
    integer :: log_limited_steps = 0
    integer :: cauchy_merit_fallback_steps = 0
    integer :: trust_region_steps = 0
    integer :: trust_region_rejections = 0
    real(dp) :: initial_residual_rms = huge(1.0_dp)
    real(dp) :: final_residual_rms = huge(1.0_dp)
    real(dp) :: final_update_rms = huge(1.0_dp)
    real(dp) :: minimum_accepted_step = 1.0_dp
    real(dp) :: maximum_armijo_violation = 0.0_dp
    real(dp) :: final_update_by_field(NDOF_NODE) = huge(1.0_dp)
    type(DiagnosticsType) :: final_diagnostics
    real(dp), allocatable :: residual_history(:)
    real(dp), allocatable :: update_history(:)
    real(dp), allocatable :: step_history(:)
  end type NewtonReportType

  public :: solve_single_load_newton
  public :: dogleg_direction

contains

  subroutine solve_single_load_newton(mesh,materials,load,context,boundary_conditions, &
                                      field_scales,options,state,report,status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(LoadStateType), intent(in) :: load
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    type(NonlinearOptionsType), intent(in) :: options
    type(StateType), intent(inout) :: state
    type(NewtonReportType), intent(out) :: report
    type(SolverStatusType), intent(out) :: status
    type(StateType) :: trial_state
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: system,trial_system
    type(DiagnosticsType) :: diagnostics,trial_diagnostics
    type(PardisoDiagnosticsType) :: linear_diagnostics
    type(SolverStatusType) :: local_status,last_trial_status
    logical, allocatable :: fixed(:)
    real(dp), allocatable :: variable_scales(:),residual_scales(:)
    real(dp), allocatable :: scaled_values(:),scaled_rhs(:),rhs(:)
    real(dp), allocatable :: scaled_increment(:),increment(:),raw_increment(:), &
      trial_scaled_residual(:)
    real(dp), allocatable :: scaled_direction(:),scaled_jacobian_action(:),scaled_residual(:)
    real(dp), allocatable :: scaled_merit_gradient(:),scaled_gradient_action(:)
    real(dp) :: residual_rms,trial_residual_rms,full_update_rms,accepted_update_rms
    real(dp) :: base_merit,trial_merit,armijo_bound,armijo_violation,step,merit_slope
    real(dp) :: component_merit_slope,raw_log_density_update,cauchy_log_factor
    real(dp) :: raw_eta_n_update,raw_eta_p_update,cauchy_alpha
    real(dp) :: trust_radius,predicted_reduction,actual_reduction,trust_ratio,direction_norm
    integer :: iteration,backtracks,equation_count,trust_strategy
    logical :: accepted,log_step_limited

    call clear_status(status)
    report=NewtonReportType()
    call validate_options(field_scales,options,status)
    if (.not. status%is_ok()) return
    equation_count=context%number_of_equations
    if (equation_count /= NDOF_NODE*mesh%number_of_nodes()) then
      call set_status(status,STATUS_INVALID_ARGUMENT, &
        'Newton sparse context does not match the mesh.')
      return
    end if
    call project_dirichlet_state(state,boundary_conditions,load,mesh%number_of_nodes(), &
      local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call build_fixed_dof_mask(boundary_conditions,mesh%number_of_nodes(),fixed,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    if (count(.not. fixed) < 1) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Newton system has no free degrees of freedom.')
      return
    end if
    call build_variable_scaling(mesh%number_of_nodes(),field_scales,variable_scales,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if

    allocate(trial_state%x(equation_count),rhs(equation_count), &
      scaled_increment(equation_count),increment(equation_count),raw_increment(equation_count), &
      trial_scaled_residual(equation_count),scaled_direction(equation_count), &
      scaled_jacobian_action(equation_count),scaled_residual(equation_count), &
      scaled_merit_gradient(equation_count),scaled_gradient_action(equation_count))
    allocate(report%residual_history(options%max_newton+1), &
      report%update_history(options%max_newton),report%step_history(options%max_newton))
    report%residual_history=-1.0_dp
    report%update_history=-1.0_dp
    report%step_history=-1.0_dp
    accepted_update_rms=huge(1.0_dp)
    request%need_jacobian=.true.
    trust_radius=options%trust_initial_radius

    do iteration=0,options%max_newton
      call assemble_global_system(mesh,materials,state,load,context,boundary_conditions, &
        request,system,diagnostics,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      rhs=-system%residual
      call equilibrate_csr_system(context,system%jacobian_values,rhs,variable_scales, &
        scaled_values,scaled_rhs,residual_scales,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      residual_rms=free_rms(scaled_rhs,fixed)
      report%residual_history(iteration+1)=residual_rms
      report%history_count=iteration+1
      report%final_residual_rms=residual_rms
      report%final_diagnostics=diagnostics
      if (options%verbose) then
        write(output_unit,'(a,i0,a,3(es12.4,a))') 'newton_iteration=',iteration, &
          ', voltage_factor=',load%voltage_factor,', joule_factor=',load%joule_factor, &
          ', residual=',residual_rms,','
        flush(output_unit)
      end if
      if (iteration == 0) report%initial_residual_rms=residual_rms

      if (residual_rms <= options%residual_tolerance .and. &
          (iteration == 0 .or. accepted_update_rms <= options%update_tolerance)) then
        report%converged=.true.
        report%iterations=iteration
        if (iteration == 0) accepted_update_rms=0.0_dp
        report%final_update_rms=accepted_update_rms
        return
      end if
      if (iteration == options%max_newton) exit

      ! Optional early rejection of a stalled LOAD ATTEMPT, never acceptance.
      ! Continuation restores the last accepted state and reduces its load step.
      ! Scaled residuals change with row equilibration, so this is a work-limit
      ! heuristic only; it makes no claim about absence of a physical solution.
      if(options%stagnation_window>0 .and.residual_rms>options%residual_tolerance) then
        if(iteration>=options%stagnation_window) then
          if(residual_rms>=0.99_dp*report%residual_history(iteration+1-options%stagnation_window)) then
            call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED, &
              'Newton residual stagnation: reject this load attempt and cut back.')
            return
          end if
        end if
      end if

      call solve_csr_pardiso(equation_count,context%ia,context%ja,scaled_values, &
        scaled_rhs,scaled_increment,linear_diagnostics,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      report%total_linear_solves=report%total_linear_solves+1
      report%maximum_perturbed_pivots=max(report%maximum_perturbed_pivots, &
        linear_diagnostics%perturbed_pivots)
      call unscale_solution(variable_scales,scaled_increment,increment,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      where(fixed) increment=0.0_dp
      raw_increment=increment
      raw_log_density_update=maximum_log_density_increment(raw_increment)
      raw_eta_n_update=maximum_field_increment(raw_increment,F_ETA_N)
      raw_eta_p_update=maximum_field_increment(raw_increment,F_ETA_P)
      if (options%trust_region) then
        ! The merit and row scales stay frozen during all trials for this Jacobian.
        ! Rejected trials reuse the same factorization, rather than another costly solve.
        base_merit=0.5_dp*residual_rms*residual_rms
        scaled_residual=system%residual/residual_scales
        where(fixed) scaled_increment=0.0_dp
        call csr_transpose_matvec_free(context,scaled_values,scaled_residual,fixed,scaled_merit_gradient)
        call csr_matvec(context,scaled_values,scaled_merit_gradient,scaled_gradient_action,local_status)
        if (.not.local_status%is_ok()) then
          call copy_failure(local_status,status);return
        end if
        cauchy_alpha=free_dot(scaled_merit_gradient,scaled_merit_gradient,fixed)/ &
          max(free_dot(scaled_gradient_action,scaled_gradient_action,fixed),tiny(1.0_dp))
        scaled_merit_gradient=-cauchy_alpha*scaled_merit_gradient
        accepted=.false.;backtracks=0;trust_strategy=0
        request%need_jacobian=.false.
        do
          if(trust_strategy==0) then
            ! Preserve the effective legacy candidate when its linear model
            ! predicts descent. Never accept clipping on faith: test its model.
            increment=raw_increment
            call limit_log_density_increment(increment,options%maximum_log_density_update,log_step_limited)
            scaled_direction=increment/variable_scales
            scaled_direction=scaled_direction*min(1.0_dp,trust_radius/ &
              max(free_rms(scaled_direction,fixed),tiny(1.0_dp)))
          else if(trust_strategy==2) then
            call dogleg_direction(scaled_increment,scaled_merit_gradient,fixed,trust_radius,scaled_direction)
          else
            ! Try the coupled Newton direction next. A Cauchy-heavy dogleg
            ! can crawl in ill-conditioned multiphysics systems despite good rho.
            scaled_direction=scaled_increment*min(1.0_dp,trust_radius/ &
              max(free_rms(scaled_increment,fixed),tiny(1.0_dp)))
          end if
          increment=variable_scales*scaled_direction
          ! One scalar cap preserves direction; component clipping need not preserve descent.
          cauchy_log_factor=min(1.0_dp,options%maximum_log_density_update/ &
            max(maximum_log_density_increment(increment),tiny(1.0_dp)))
          increment=cauchy_log_factor*increment
          scaled_direction=increment/variable_scales
          direction_norm=free_rms(scaled_direction,fixed)
          call csr_matvec(context,scaled_values,scaled_direction,scaled_jacobian_action,local_status)
          if (.not.local_status%is_ok()) then
            call copy_failure(local_status,status);return
          end if
          predicted_reduction=-free_dot(scaled_residual,scaled_jacobian_action,fixed)- &
            0.5_dp*free_dot(scaled_jacobian_action,scaled_jacobian_action,fixed)
          trust_ratio=-huge(1.0_dp)
          if (ieee_is_finite(predicted_reduction).and.predicted_reduction>0.0_dp) then
            trial_state%x=state%x+increment
            call project_dirichlet_state(trial_state,boundary_conditions,load,mesh%number_of_nodes(),local_status)
            if(local_status%is_ok()) call assemble_global_system(mesh,materials,trial_state,load,context, &
              boundary_conditions,request,trial_system,trial_diagnostics,local_status)
            if(local_status%is_ok()) then
              trial_scaled_residual=trial_system%residual/residual_scales
              trial_residual_rms=free_rms(trial_scaled_residual,fixed)
              trial_merit=0.5_dp*trial_residual_rms*trial_residual_rms
              actual_reduction=base_merit-trial_merit
              if(ieee_is_finite(actual_reduction)) trust_ratio=actual_reduction/predicted_reduction
            end if
          end if
          accepted=ieee_is_finite(trust_ratio).and.trust_ratio>1.0e-4_dp
          if(options%verbose) then
            write(output_unit,'(a,l1,a,i0,a,4(es12.4,a))') 'newton_trust_trial, accepted=',accepted, &
              ', strategy=',trust_strategy, &
              ', radius=',trust_radius,', step_norm=',direction_norm,', predicted=',predicted_reduction, &
              ', rho=',trust_ratio,','
            flush(output_unit)
          end if
          if(.not.accepted) then
            backtracks=backtracks+1
            report%trust_region_rejections=report%trust_region_rejections+1
            if(trust_strategy<2) then
              trust_strategy=trust_strategy+1
              cycle
            end if
          end if
          if(trust_ratio<0.25_dp.or..not.ieee_is_finite(trust_ratio)) then
            trust_radius=0.25_dp*min(trust_radius,max(direction_norm,options%trust_minimum_radius))
          else if(trust_ratio>0.75_dp.and.direction_norm>=0.8_dp*trust_radius) then
            trust_radius=min(2.0_dp*trust_radius,options%trust_maximum_radius)
          end if
          if(accepted) exit
          trust_strategy=0
          if(trust_radius<options%trust_minimum_radius.or.backtracks>=40) then
            call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED, &
              'Trust region exhausted without true residual reduction: reject load attempt.')
            return
          end if
        end do
        request%need_jacobian=.true.
        state%x=trial_state%x
        accepted_update_rms=direction_norm
        call scaled_update_by_field(increment,variable_scales,report%final_update_by_field)
        report%iterations=iteration+1
        report%trust_region_steps=report%trust_region_steps+1
        report%total_backtracks=report%total_backtracks+backtracks
        if(cauchy_log_factor<1.0_dp) report%log_limited_steps=report%log_limited_steps+1
        report%update_history(iteration+1)=accepted_update_rms
        report%step_history(iteration+1)=1.0_dp
        report%final_update_rms=accepted_update_rms
        cycle
      end if
      if (options%global_log_step_limit) then
        log_step_limited=raw_log_density_update>options%maximum_log_density_update
        if(log_step_limited) increment=increment*(options%maximum_log_density_update/raw_log_density_update)
        if(options%verbose .and.log_step_limited) then
          write(output_unit,'(a,3(es12.4,a))') 'newton_global_log_limit, raw_eta_n=',raw_eta_n_update, &
            ', raw_eta_p=',raw_eta_p_update,', scale=',options%maximum_log_density_update/raw_log_density_update,','
          flush(output_unit)
        end if
      else
        call limit_log_density_increment(increment,options%maximum_log_density_update,log_step_limited)
      end if
      if (log_step_limited) report%log_limited_steps=report%log_limited_steps+1

      base_merit=0.5_dp*residual_rms*residual_rms
      scaled_direction=increment/variable_scales
      scaled_residual=system%residual/residual_scales
      call csr_matvec(context,scaled_values,scaled_direction,scaled_jacobian_action,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      merit_slope=free_dot(scaled_residual,scaled_jacobian_action,fixed)
      component_merit_slope=merit_slope
      cauchy_log_factor=1.0_dp
      cauchy_alpha=0.0_dp
      if (log_step_limited .and. &
          (.not. ieee_is_finite(merit_slope) .or. merit_slope >= 0.0_dp)) then
        ! Independent clipping of eta components can rotate an otherwise valid Newton
        ! direction out of the descent half-space.  Fall back to the Cauchy direction
        ! for the same frozen scaled-residual merit.  Unlike globally shrinking the
        ! raw Newton step, this direction does not inherit an ill-conditioned inverse
        ! update from a nearly inactive minority-carrier block.
        call csr_transpose_matvec_free(context,scaled_values,scaled_residual,fixed, &
          scaled_merit_gradient)
        call csr_matvec(context,scaled_values,scaled_merit_gradient, &
          scaled_gradient_action,local_status)
        if (.not. local_status%is_ok()) then
          call copy_failure(local_status,status); return
        end if
        cauchy_alpha=free_dot(scaled_merit_gradient,scaled_merit_gradient,fixed)/ &
          max(free_dot(scaled_gradient_action,scaled_gradient_action,fixed),tiny(1.0_dp))
        scaled_direction=-cauchy_alpha*scaled_merit_gradient
        increment=variable_scales*scaled_direction
        where(fixed) increment=0.0_dp
        cauchy_log_factor=min(1.0_dp,options%maximum_log_density_update/ &
          max(maximum_log_density_increment(increment),tiny(1.0_dp)))
        increment=cauchy_log_factor*increment
        scaled_direction=increment/variable_scales
        call csr_matvec(context,scaled_values,scaled_direction,scaled_jacobian_action, &
          local_status)
        if (.not. local_status%is_ok()) then
          call copy_failure(local_status,status); return
        end if
        merit_slope=free_dot(scaled_residual,scaled_jacobian_action,fixed)
        report%cauchy_merit_fallback_steps=report%cauchy_merit_fallback_steps+1
        if (options%verbose) then
          write(output_unit,'(a,7(es12.4,a))') 'newton_descent_fallback=CAUCHY_MERIT, raw_eta_n=', &
            raw_eta_n_update,', raw_eta_p=',raw_eta_p_update, &
            ', component_slope=',component_merit_slope,', alpha=',cauchy_alpha, &
            ', log_scale=',cauchy_log_factor,', fallback_slope=',merit_slope, &
            ', raw_max=',raw_log_density_update,','
          flush(output_unit)
        end if
      end if
      if (.not. ieee_is_finite(merit_slope) .or. merit_slope >= 0.0_dp) then
        call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED, &
          'Safeguarded Newton direction is not a descent direction for the true residual merit.')
        return
      end if
      full_update_rms=scaled_free_rms(increment,variable_scales,fixed)
      step=1.0_dp
      backtracks=0
      accepted=.false.
      call clear_status(last_trial_status)
      request%need_jacobian=.false.
      do
        trial_state%x=state%x+step*increment
        call project_dirichlet_state(trial_state,boundary_conditions,load, &
          mesh%number_of_nodes(),local_status)
        if (local_status%is_ok()) then
          call assemble_global_system(mesh,materials,trial_state,load,context, &
            boundary_conditions,request,trial_system,trial_diagnostics,local_status)
        end if
        if (local_status%is_ok()) then
          trial_scaled_residual=-trial_system%residual/residual_scales
          trial_residual_rms=free_rms(trial_scaled_residual,fixed)
          trial_merit=0.5_dp*trial_residual_rms*trial_residual_rms
          armijo_bound=base_merit+options%armijo_c1*step*merit_slope
          if (ieee_is_finite(trial_merit) .and. trial_merit <= armijo_bound) then
            armijo_violation=max(0.0_dp,(trial_merit-armijo_bound)/max(base_merit,tiny(1.0_dp)))
            report%maximum_armijo_violation=max(report%maximum_armijo_violation, &
              armijo_violation)
            accepted=.true.
            exit
          end if
        else
          last_trial_status=local_status
        end if
        step=step*options%backtrack_factor
        backtracks=backtracks+1
        if (step < options%minimum_step) exit
      end do
      request%need_jacobian=.true.

      if (.not. accepted) then
        if (.not. last_trial_status%is_ok()) then
          call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED, &
            'Armijo line search reached minimum step; last trial failed: '// &
            trim(last_trial_status%message))
        else
          call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED, &
            'Armijo line search reached the frozen minimum step without residual decrease.')
        end if
        return
      end if

      state%x=trial_state%x
      accepted_update_rms=step*full_update_rms
      ! Avoid a large compiler-generated stack temporary for step*increment.
      call scaled_update_by_field(increment,variable_scales, &
        report%final_update_by_field)
      report%final_update_by_field=step*report%final_update_by_field
      report%iterations=iteration+1
      report%total_backtracks=report%total_backtracks+backtracks
      report%minimum_accepted_step=min(report%minimum_accepted_step,step)
      report%update_history(iteration+1)=accepted_update_rms
      report%step_history(iteration+1)=step
      report%final_update_rms=accepted_update_rms
    end do

    call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED, &
      'Newton iteration reached the frozen maximum without satisfying both tolerances.')
  end subroutine solve_single_load_newton

  subroutine validate_options(field_scales,options,status)
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    type(NonlinearOptionsType), intent(in) :: options
    type(SolverStatusType), intent(out) :: status

    call clear_status(status)
    if (any(field_scales <= 0.0_dp) .or. .not. all(ieee_is_finite(field_scales))) then
      call set_status(status,STATUS_INVALID_ARGUMENT, &
        'Newton field scales must be finite and positive.')
    else if (options%max_newton < 1 .or. options%stagnation_window < 0 .or. &
             options%stagnation_window > options%max_newton .or. options%residual_tolerance <= 0.0_dp .or. &
             options%update_tolerance <= 0.0_dp) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Newton iteration limits are invalid.')
    else if (options%armijo_c1 <= 0.0_dp .or. options%armijo_c1 >= 0.5_dp .or. &
             options%backtrack_factor <= 0.0_dp .or. options%backtrack_factor >= 1.0_dp .or. &
             options%minimum_step <= 0.0_dp .or. options%minimum_step > 1.0_dp .or. &
             options%maximum_log_density_update <= 0.0_dp .or. &
             .not. ieee_is_finite(options%maximum_log_density_update)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Newton Armijo options are invalid.')
    end if
    if (.not.status%is_ok()) return
    if (options%trust_region) then
      if (.not.all(ieee_is_finite([options%trust_initial_radius,options%trust_maximum_radius, &
          options%trust_minimum_radius])).or.options%trust_minimum_radius<=0.0_dp.or. &
          options%trust_initial_radius<options%trust_minimum_radius.or. &
          options%trust_maximum_radius<options%trust_initial_radius) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'Newton trust region radii are invalid.')
      end if
    end if
  end subroutine validate_options

  subroutine dogleg_direction(newton,cauchy,fixed,radius,direction)
    real(dp), intent(in) :: newton(:),cauchy(:),radius
    logical, intent(in) :: fixed(:)
    real(dp), intent(out) :: direction(:)
    real(dp) :: norm_n,norm_c,a,b,c,tau,discriminant
    norm_n=free_rms(newton,fixed);norm_c=free_rms(cauchy,fixed)
    if(norm_n<=radius) then
      direction=newton
    else if(norm_c>=radius) then
      direction=cauchy*(radius/max(norm_c,tiny(1.0_dp)))
    else
      direction=newton-cauchy
      a=free_dot(direction,direction,fixed)
      b=free_dot(cauchy,direction,fixed)
      c=norm_c**2-radius**2
      discriminant=sqrt(max(0.0_dp,b*b-a*c))
      if(b>=0.0_dp) then
        tau=-c/max(b+discriminant,tiny(1.0_dp))
      else
        tau=(-b+discriminant)/max(a,tiny(1.0_dp))
      end if
      direction=cauchy+max(0.0_dp,min(1.0_dp,tau))*direction
    end if
    where(fixed) direction=0.0_dp
  end subroutine dogleg_direction

  real(dp) function free_rms(vector,fixed) result(value)
    real(dp), intent(in) :: vector(:)
    logical, intent(in) :: fixed(:)
    real(dp) :: square_sum
    integer :: equation,free_count

    square_sum=0.0_dp
    free_count=0
    do equation=1,size(vector)
      if (fixed(equation)) cycle
      square_sum=square_sum+vector(equation)*vector(equation)
      free_count=free_count+1
    end do
    value=sqrt(square_sum/real(free_count,dp))
  end function free_rms

  real(dp) function scaled_free_rms(vector,scales,fixed) result(value)
    real(dp), intent(in) :: vector(:),scales(:)
    logical, intent(in) :: fixed(:)
    real(dp) :: square_sum
    integer :: equation,free_count

    square_sum=0.0_dp
    free_count=0
    do equation=1,size(vector)
      if (fixed(equation)) cycle
      square_sum=square_sum+(vector(equation)/scales(equation))**2
      free_count=free_count+1
    end do
    value=sqrt(square_sum/real(free_count,dp))
  end function scaled_free_rms

  subroutine scaled_update_by_field(vector,scales,field_rms)
    real(dp), intent(in) :: vector(:),scales(:)
    real(dp), intent(out) :: field_rms(NDOF_NODE)
    integer :: equation,field,field_count(NDOF_NODE)

    field_rms=0.0_dp
    field_count=0
    do equation=1,size(vector)
      field=mod(equation-1,NDOF_NODE)+1
      field_rms(field)=field_rms(field)+(vector(equation)/scales(equation))**2
      field_count(field)=field_count(field)+1
    end do
    do field=1,NDOF_NODE
      field_rms(field)=sqrt(field_rms(field)/real(field_count(field),dp))
    end do
  end subroutine scaled_update_by_field

  subroutine limit_log_density_increment(increment,maximum_change,limited)
    real(dp), intent(inout) :: increment(:)
    real(dp), intent(in) :: maximum_change
    logical, intent(out) :: limited
    integer :: equation,field

    limited=.false.
    do equation=1,size(increment)
      field=mod(equation-1,NDOF_NODE)+1
      if (field == F_ETA_N .or. field == F_ETA_P) then
        if (abs(increment(equation)) > maximum_change) then
          increment(equation)=sign(maximum_change,increment(equation))
          limited=.true.
        end if
      end if
    end do
  end subroutine limit_log_density_increment

  real(dp) function maximum_log_density_increment(increment) result(maximum_change)
    real(dp), intent(in) :: increment(:)
    integer :: equation,field

    maximum_change=0.0_dp
    do equation=1,size(increment)
      field=mod(equation-1,NDOF_NODE)+1
      if (field == F_ETA_N .or. field == F_ETA_P) then
        maximum_change=max(maximum_change,abs(increment(equation)))
      end if
    end do
  end function maximum_log_density_increment

  real(dp) function maximum_field_increment(increment,target_field) result(maximum_change)
    real(dp), intent(in) :: increment(:)
    integer, intent(in) :: target_field
    integer :: equation,field

    maximum_change=0.0_dp
    do equation=1,size(increment)
      field=mod(equation-1,NDOF_NODE)+1
      if (field == target_field) maximum_change=max(maximum_change,abs(increment(equation)))
    end do
  end function maximum_field_increment

  subroutine csr_transpose_matvec_free(context,values,vector,fixed,result)
    type(SparseAssemblyContextType), intent(in) :: context
    real(dp), intent(in) :: values(:),vector(:)
    logical, intent(in) :: fixed(:)
    real(dp), intent(out) :: result(:)
    integer :: row,position,column

    result=0.0_dp
    do row=1,context%number_of_equations
      if (fixed(row)) cycle
      do position=context%ia(row),context%ia(row+1)-1
        column=context%ja(position)
        if (.not. fixed(column)) result(column)=result(column)+values(position)*vector(row)
      end do
    end do
  end subroutine csr_transpose_matvec_free

  real(dp) function free_dot(vector_a,vector_b,fixed) result(value)
    real(dp), intent(in) :: vector_a(:),vector_b(:)
    logical, intent(in) :: fixed(:)
    real(dp) :: product_sum
    integer :: equation,free_count

    product_sum=0.0_dp
    free_count=0
    do equation=1,size(vector_a)
      if (fixed(equation)) cycle
      product_sum=product_sum+vector_a(equation)*vector_b(equation)
      free_count=free_count+1
    end do
    value=product_sum/real(free_count,dp)
  end function free_dot

  subroutine copy_failure(source,destination)
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,source%code,trim(source%message))
  end subroutine copy_failure

end module newton_solver_m

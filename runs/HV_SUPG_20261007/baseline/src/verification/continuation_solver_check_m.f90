module continuation_solver_check_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE
  use input_schema_m, only : CaseConfigType, MODEL_M2, PROBLEM_REFERENCE_BAR, &
    MATERIAL_TPS_G08_REFERENCE
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, ContinuationOptionsType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use scaling_m, only : build_variable_scaling
  use reference_bar_problem_m, only : initialize_reference_bar_problem
  use continuation_solver_m, only : ContinuationReportType, solve_load_continuation
  use global_assembly_check_m, only : build_two_element_mesh
  implicit none
  private

  public :: verify_continuation_solver_g08

contains

  subroutine verify_continuation_solver_g08(status)
    type(SolverStatusType), intent(out) :: status
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    type(MaterialType), allocatable :: materials(:)
    type(StateType) :: full_state,cutback_state,restart_state
    type(LoadStateType) :: full_load,cutback_load,restart_load,target_load
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(SparseAssemblyContextType) :: context
    type(ContinuationReportType) :: full_report,cutback_report,restart_report
    type(ContinuationOptionsType) :: cutback_options
    type(SolverStatusType) :: local_status
    real(dp) :: restart_difference,cutback_difference,maximum_residual,maximum_update
    logical :: passed

    call clear_status(status)
    call define_g08_config(config)
    call build_two_element_mesh(mesh)
    mesh%coordinates=1.0e-6_dp*mesh%coordinates
    call initialize_reference_bar_problem(config,mesh,materials,full_state,full_load, &
      boundary_conditions,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G08 problem setup failed: ',local_status,status); return
    end if
    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G08 sparse setup failed: ',local_status,status); return
    end if
    target_load%voltage_factor=config%target_voltage_factor
    target_load%joule_factor=config%target_joule_factor
    call solve_load_continuation(mesh,materials,context,boundary_conditions, &
      config%field_scales,config%nonlinear_options,config%continuation_options, &
      target_load,full_state,full_load,full_report,local_status)
    if (.not. local_status%is_ok()) then
      write(*,'(a)') 'verification_id=G08_continuation_production'
      call print_report('full',full_report)
      call copy_failure('G08 full continuation failed: ',local_status,status); return
    end if
    if (.not. allocated(full_report%transition_state%x)) then
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G08 voltage-complete restart checkpoint is missing.')
      return
    end if

    allocate(cutback_state%x(size(full_state%x)))
    cutback_state%x=0.0_dp
    cutback_load%voltage_factor=0.0_dp
    cutback_load%joule_factor=0.0_dp
    cutback_options=config%continuation_options
    cutback_options%maximum_newton_per_attempt=7
    call solve_load_continuation(mesh,materials,context,boundary_conditions, &
      config%field_scales,config%nonlinear_options,cutback_options,target_load, &
      cutback_state,cutback_load,cutback_report,local_status)
    if (.not. local_status%is_ok()) then
      write(*,'(a)') 'verification_id=G08_continuation_production'
      call print_report('full',full_report)
      call print_report('cutback',cutback_report)
      call copy_failure('G08 cutback stress path failed: ',local_status,status); return
    end if

    allocate(restart_state%x(size(full_report%transition_state%x)))
    restart_state%x=full_report%transition_state%x
    restart_load=full_report%transition_load
    call solve_load_continuation(mesh,materials,context,boundary_conditions, &
      config%field_scales,config%nonlinear_options,config%continuation_options, &
      target_load,restart_state,restart_load,restart_report,local_status)
    if (.not. local_status%is_ok()) then
      write(*,'(a)') 'verification_id=G08_continuation_production'
      call print_report('full',full_report)
      call print_report('restart',restart_report)
      call copy_failure('G08 restart continuation failed: ',local_status,status); return
    end if

    restart_difference=scaled_state_difference(mesh,config%field_scales, &
      full_state,restart_state,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G08 restart comparison failed: ',local_status,status); return
    end if
    cutback_difference=scaled_state_difference(mesh,config%field_scales, &
      full_state,cutback_state,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G08 cutback comparison failed: ',local_status,status); return
    end if
    maximum_residual=max(maxval(full_report%residual_rms(1:full_report%accepted_points)), &
      maxval(cutback_report%residual_rms(1:cutback_report%accepted_points)), &
      maxval(restart_report%residual_rms(1:restart_report%accepted_points)))
    maximum_update=max(maxval(full_report%update_rms(1:full_report%accepted_points)), &
      maxval(cutback_report%update_rms(1:cutback_report%accepted_points)), &
      maxval(restart_report%update_rms(1:restart_report%accepted_points)))
    passed=full_report%converged .and. restart_report%converged .and. &
      cutback_report%converged .and. &
      full_load%voltage_factor == target_load%voltage_factor .and. &
      full_load%joule_factor == target_load%joule_factor .and. &
      restart_load%voltage_factor == target_load%voltage_factor .and. &
      restart_load%joule_factor == target_load%joule_factor .and. &
      cutback_load%voltage_factor == target_load%voltage_factor .and. &
      cutback_load%joule_factor == target_load%joule_factor .and. &
      cutback_report%rejected_attempts >= 1 .and. &
      cutback_report%maximum_checkpoint_change_on_rejection == 0.0_dp .and. &
      maximum_residual < config%nonlinear_options%residual_tolerance .and. &
      maximum_update < config%nonlinear_options%update_tolerance .and. &
      max(full_report%maximum_balance_defect,cutback_report%maximum_balance_defect, &
        restart_report%maximum_balance_defect) < 1.0e-8_dp .and. &
      max(full_report%maximum_boundary_error,cutback_report%maximum_boundary_error, &
        restart_report%maximum_boundary_error) < 1.0e-14_dp .and. &
      report_is_physically_admissible(full_report) .and. &
      report_is_physically_admissible(cutback_report) .and. &
      report_is_physically_admissible(restart_report) .and. &
      restart_difference < 1.0e-12_dp .and. cutback_difference < 1.0e-12_dp

    write(*,'(a)') 'verification_id=G08_continuation_production'
    call print_report('full',full_report)
    call print_report('cutback',cutback_report)
    call print_report('restart',restart_report)
    call print_metric('g08_restart_scaled_solution_difference',restart_difference)
    call print_metric('g08_cutback_scaled_solution_difference',cutback_difference)
    call print_metric('g08_maximum_accepted_residual_rms',maximum_residual)
    call print_metric('g08_maximum_accepted_update_rms',maximum_update)
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G08 path, cutback, restart, balance, or admissibility requirement failed.')
    end if
  end subroutine verify_continuation_solver_g08

  subroutine define_g08_config(config)
    type(CaseConfigType), intent(out) :: config

    config=CaseConfigType()
    config%schema_version=2
    config%case_id='G08_TWO_ELEMENT_REFERENCE_BAR'
    config%model_level=MODEL_M2
    config%mesh_file='IN_MEMORY'
    config%problem_type=PROBLEM_REFERENCE_BAR
    config%material_profile=MATERIAL_TPS_G08_REFERENCE
    config%target_voltage=-10.0_dp
    config%target_voltage_factor=1.0_dp
    config%target_joule_factor=1.0_dp
    config%field_scales=[1.0e-10_dp,1.0e-10_dp,1.0e-10_dp,10.0_dp, &
      1.0e-1_dp,1.0e-1_dp,10.0_dp]
    config%nonlinear_options%max_newton=25
    config%nonlinear_options%residual_tolerance=1.0e-9_dp
    config%nonlinear_options%update_tolerance=1.0e-9_dp
    config%nonlinear_options%armijo_c1=1.0e-4_dp
    config%nonlinear_options%backtrack_factor=0.5_dp
    config%nonlinear_options%minimum_step=1.0e-8_dp
    config%nonlinear_options%maximum_log_density_update=1.0_dp
    config%continuation_options%initial_step=1.0_dp
    config%continuation_options%minimum_step=0.0625_dp
    config%continuation_options%maximum_step=1.0_dp
    config%continuation_options%cutback_factor=0.5_dp
    config%continuation_options%growth_factor=2.0_dp
    config%continuation_options%maximum_attempts=64
    config%continuation_options%maximum_newton_per_attempt=8
    config%continuation_options%fast_iteration_threshold=3
    config%continuation_options%slow_iteration_threshold=7
  end subroutine define_g08_config

  logical function report_is_physically_admissible(report) result(admissible)
    type(ContinuationReportType), intent(in) :: report
    integer :: point

    admissible=.true.
    do point=1,report%accepted_points
      admissible=admissible .and. report%minimum_temperature(point) > 0.0_dp .and. &
        report%minimum_n(point) > 0.0_dp .and. report%minimum_p(point) > 0.0_dp .and. &
        ieee_is_finite(report%minimum_temperature(point)) .and. &
        ieee_is_finite(report%maximum_temperature(point)) .and. &
        ieee_is_finite(report%minimum_n(point)) .and. ieee_is_finite(report%maximum_n(point)) .and. &
        ieee_is_finite(report%minimum_p(point)) .and. ieee_is_finite(report%maximum_p(point))
    end do
  end function report_is_physically_admissible

  real(dp) function scaled_state_difference(mesh,field_scales,state_a,state_b,status) &
    result(value)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    type(StateType), intent(in) :: state_a,state_b
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: variable_scales(:)
    integer :: equation

    call build_variable_scaling(mesh%number_of_nodes(),field_scales,variable_scales,status)
    if (.not. status%is_ok()) then
      value=huge(1.0_dp); return
    end if
    value=0.0_dp
    do equation=1,size(state_a%x)
      value=value+((state_a%x(equation)-state_b%x(equation))/variable_scales(equation))**2
    end do
    value=sqrt(value/real(size(state_a%x),dp))
  end function scaled_state_difference

  subroutine print_report(label,report)
    character(len=*), intent(in) :: label
    type(ContinuationReportType), intent(in) :: report
    integer :: point

    write(*,'(a,l1)') 'g08_'//trim(label)//'_converged=',report%converged
    write(*,'(a,i0)') 'g08_'//trim(label)//'_accepted_points=',report%accepted_points
    write(*,'(a,i0)') 'g08_'//trim(label)//'_rejected_attempts=',report%rejected_attempts
    write(*,'(a,i0)') 'g08_'//trim(label)//'_total_attempts=',report%total_attempts
    write(*,'(a,i0)') 'g08_'//trim(label)//'_total_newton_iterations=', &
      report%total_newton_iterations
    call print_metric('g08_'//trim(label)//'_minimum_accepted_step', &
      report%minimum_accepted_step)
    call print_metric('g08_'//trim(label)//'_checkpoint_change_on_rejection', &
      report%maximum_checkpoint_change_on_rejection)
    call print_metric('g08_'//trim(label)//'_maximum_balance_defect', &
      report%maximum_balance_defect)
    call print_metric('g08_'//trim(label)//'_maximum_boundary_error', &
      report%maximum_boundary_error)
    do point=1,report%accepted_points
      call print_point(label,point,report)
    end do
  end subroutine print_report

  subroutine print_point(label,point,report)
    character(len=*), intent(in) :: label
    integer, intent(in) :: point
    type(ContinuationReportType), intent(in) :: report
    character(len=16) :: point_text

    write(point_text,'(i0)') point
    write(*,'(a,i0)') 'g08_'//trim(label)//'_phase_'//trim(point_text)//'=', &
      report%phase(point)
    write(*,'(a,i0)') 'g08_'//trim(label)//'_newton_iterations_'//trim(point_text)//'=', &
      report%newton_iterations(point)
    call print_metric('g08_'//trim(label)//'_voltage_factor_'//trim(point_text), &
      report%voltage_factor(point))
    call print_metric('g08_'//trim(label)//'_joule_factor_'//trim(point_text), &
      report%joule_factor(point))
    call print_metric('g08_'//trim(label)//'_residual_'//trim(point_text), &
      report%residual_rms(point))
    call print_metric('g08_'//trim(label)//'_update_'//trim(point_text), &
      report%update_rms(point))
    call print_metric('g08_'//trim(label)//'_balance_'//trim(point_text), &
      maxval(report%balance(:,point)))
  end subroutine print_point

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

  subroutine copy_failure(prefix,source,destination)
    character(len=*), intent(in) :: prefix
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,STATUS_VERIFICATION_FAILED,trim(prefix)//trim(source%message))
  end subroutine copy_failure

end module continuation_solver_check_m

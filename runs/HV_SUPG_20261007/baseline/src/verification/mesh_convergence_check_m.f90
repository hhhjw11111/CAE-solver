module mesh_convergence_check_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use input_schema_m, only : CaseConfigType, PROBLEM_CENTER_HOLE_PLATE
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_VERIFICATION_FAILED
  use native_input_m, only : read_native_case
  use mesh_checks_m, only : validate_mesh
  use boundary_topology_m, only : BoundaryTopologyAuditType, audit_boundary_topology
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use reference_bar_problem_m, only : initialize_voltage_driven_problem
  use continuation_solver_m, only : ContinuationReportType, solve_load_continuation
  use field_metrics_m, only : HoleFieldMetricsType, evaluate_hole_field_metrics
  implicit none
  private

  integer, parameter :: LEVEL_COUNT=3,METRIC_COUNT=6
  real(dp), parameter :: ANNULUS_WIDTH=30.0e-9_dp
  integer, parameter :: EXPECTED_NODES(LEVEL_COUNT)=[304,624,1056]
  integer, parameter :: EXPECTED_ELEMENTS(LEVEL_COUNT)=[32,72,128]
  integer, parameter :: EXPECTED_THETA(LEVEL_COUNT)=[16,24,32]
  integer, parameter :: EXPECTED_CONTACT_NODES(LEVEL_COUNT)=[23,33,43]

  type :: G10CaseRunType
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    type(MaterialType), allocatable :: materials(:)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(SparseAssemblyContextType) :: context
    type(ContinuationReportType) :: report
    type(BoundaryTopologyAuditType) :: topology
    type(HoleFieldMetricsType) :: metrics
    real(dp) :: minimum_determinant=0.0_dp
    real(dp) :: volume=0.0_dp
    real(dp) :: volume_error=huge(1.0_dp)
    integer :: hole_dirichlet_overlap=0
  end type G10CaseRunType

  public :: verify_mesh_convergence_g10

contains

  subroutine verify_mesh_convergence_g10(case_paths,status)
    character(len=*), intent(in) :: case_paths(LEVEL_COUNT)
    type(SolverStatusType), intent(out) :: status
    type(G10CaseRunType) :: runs(LEVEL_COUNT)
    type(SolverStatusType) :: local_status
    real(dp) :: values(METRIC_COUNT,LEVEL_COUNT),change_01(METRIC_COUNT)
    real(dp) :: change_12(METRIC_COUNT),exact_volume
    logical :: topology_pass,solver_pass,threshold_pass
    integer :: level,metric,improved_count

    call clear_status(status)
    do level=1,LEVEL_COUNT
      call prepare_solve_and_measure(case_paths(level),runs(level),local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure('G10 L'//integer_to_string(level-1)//' failed: ', &
          local_status,status)
        return
      end if
      if (trim(runs(level)%config%problem_type) /= PROBLEM_CENTER_HOLE_PLATE) then
        call set_status(status,STATUS_VERIFICATION_FAILED, &
          'G10 requires three CENTER_HOLE_PLATE cases.')
        return
      end if
      exact_volume=((runs(level)%topology%bounds_max(1)- &
        runs(level)%topology%bounds_min(1))* &
        (runs(level)%topology%bounds_max(2)-runs(level)%topology%bounds_min(2))- &
        acos(-1.0_dp)*runs(level)%config%hole_radius**2)* &
        (runs(level)%topology%bounds_max(3)-runs(level)%topology%bounds_min(3))
      runs(level)%volume_error=abs(runs(level)%volume-exact_volume)/exact_volume
      runs(level)%hole_dirichlet_overlap=count_hole_dirichlet_overlap(runs(level))
      values(:,level)=metric_values(runs(level)%metrics)
    end do

    change_01=relative_change(values(:,1),values(:,2))
    change_12=relative_change(values(:,2),values(:,3))
    improved_count=0
    do metric=1,METRIC_COUNT
      if (change_12(metric) <= change_01(metric)) improved_count=improved_count+1
    end do
    topology_pass=.true.; solver_pass=.true.
    do level=1,LEVEL_COUNT
      topology_pass=topology_pass .and. topology_passes(runs(level),level)
      solver_pass=solver_pass .and. continuation_passes(runs(level))
    end do
    threshold_pass=change_12(1) < 0.01_dp .and. change_12(2) < 0.05_dp .and. &
      change_12(3) < 0.10_dp .and. all(change_12(4:6) < 0.05_dp) .and. &
      improved_count >= 4

    write(*,'(a)') 'verification_id=G10_mesh_convergence'
    call print_metric('g10_annulus_width',ANNULUS_WIDTH)
    do level=1,LEVEL_COUNT
      call print_level_metrics(level-1,runs(level),values(:,level))
    end do
    do metric=1,METRIC_COUNT
      call print_metric('g10_'//trim(metric_name(metric))//'_change_l0_l1',change_01(metric))
      call print_metric('g10_'//trim(metric_name(metric))//'_change_l1_l2',change_12(metric))
    end do
    write(*,'(a,i0)') 'g10_improved_metric_count=',improved_count
    if (topology_pass .and. solver_pass .and. threshold_pass) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G10 topology, nonlinear, conservation, or frozen mesh-convergence requirement failed.')
    end if
  end subroutine verify_mesh_convergence_g10

  subroutine prepare_solve_and_measure(case_path,run,status)
    character(len=*), intent(in) :: case_path
    type(G10CaseRunType), intent(out) :: run
    type(SolverStatusType), intent(out) :: status
    type(LoadStateType) :: target_load

    call read_native_case(case_path,run%config,run%mesh,status)
    if (.not. status%is_ok()) return
    call validate_mesh(run%mesh,run%minimum_determinant,run%volume,status)
    if (.not. status%is_ok()) return
    call audit_boundary_topology(run%mesh,run%config%hole_radius,run%topology,status)
    if (.not. status%is_ok()) return
    call initialize_voltage_driven_problem(run%config,run%mesh,run%materials,run%state, &
      run%load,run%boundary_conditions,status)
    if (.not. status%is_ok()) return
    call initialize_sparse_context(run%mesh,run%context,status)
    if (.not. status%is_ok()) return
    target_load%voltage_factor=run%config%target_voltage_factor
    target_load%joule_factor=run%config%target_joule_factor
    call solve_load_continuation(run%mesh,run%materials,run%context, &
      run%boundary_conditions,run%config%field_scales,run%config%nonlinear_options, &
      run%config%continuation_options,target_load,run%state,run%load,run%report,status)
    if (.not. status%is_ok()) return
    call evaluate_hole_field_metrics(run%mesh,run%materials,run%state,run%load, &
      run%config%hole_radius,ANNULUS_WIDTH,run%metrics,status)
  end subroutine prepare_solve_and_measure

  logical function topology_passes(run,level) result(passed)
    type(G10CaseRunType), intent(in) :: run
    integer, intent(in) :: level
    integer :: side_faces

    side_faces=EXPECTED_THETA(level)/4
    passed=run%mesh%number_of_nodes() == EXPECTED_NODES(level) .and. &
      run%mesh%number_of_elements() == EXPECTED_ELEMENTS(level) .and. &
      run%topology%left_face_count == side_faces .and. &
      run%topology%right_face_count == side_faces .and. &
      run%topology%front_face_count == side_faces .and. &
      run%topology%back_face_count == side_faces .and. &
      run%topology%bottom_face_count == EXPECTED_ELEMENTS(level) .and. &
      run%topology%top_face_count == EXPECTED_ELEMENTS(level) .and. &
      run%topology%hole_face_count == EXPECTED_THETA(level) .and. &
      run%topology%hole_node_count == 5*EXPECTED_THETA(level) .and. &
      run%topology%unknown_face_count == 0 .and. &
      run%topology%left_node_count == EXPECTED_CONTACT_NODES(level) .and. &
      run%topology%right_node_count == EXPECTED_CONTACT_NODES(level) .and. &
      run%topology%maximum_hole_radius_error <= &
        run%topology%classification_tolerance .and. &
      run%hole_dirichlet_overlap == 0 .and. run%minimum_determinant > 0.0_dp .and. &
      run%volume_error < 1.0e-3_dp
  end function topology_passes

  logical function continuation_passes(run) result(passed)
    type(G10CaseRunType), intent(in) :: run
    integer :: point

    passed=run%report%converged .and. run%load%voltage_factor == 1.0_dp .and. &
      run%load%joule_factor == 1.0_dp .and. &
      run%report%maximum_balance_defect < 1.0e-8_dp .and. &
      run%report%maximum_boundary_error < 1.0e-14_dp
    do point=1,run%report%accepted_points
      passed=passed .and. run%report%residual_rms(point) < 1.0e-9_dp .and. &
        run%report%update_rms(point) < 1.0e-9_dp .and. &
        run%report%minimum_temperature(point) > 0.0_dp .and. &
        run%report%minimum_n(point) > 0.0_dp .and. run%report%minimum_p(point) > 0.0_dp .and. &
        ieee_is_finite(run%report%maximum_temperature(point)) .and. &
        ieee_is_finite(run%report%maximum_n(point)) .and. &
        ieee_is_finite(run%report%maximum_p(point))
    end do
  end function continuation_passes

  integer function count_hole_dirichlet_overlap(run) result(overlap)
    type(G10CaseRunType), intent(in) :: run
    logical, allocatable :: constrained(:),hole_node(:)
    real(dp) :: radius,tolerance
    integer :: boundary,node

    allocate(constrained(run%mesh%number_of_nodes()),hole_node(run%mesh%number_of_nodes()))
    constrained=.false.; hole_node=.false.
    do boundary=1,size(run%boundary_conditions)
      constrained(run%boundary_conditions(boundary)%node_ids)=.true.
    end do
    tolerance=run%topology%classification_tolerance
    do node=1,run%mesh%number_of_nodes()
      radius=hypot(run%mesh%coordinates(1,node)-run%topology%center(1), &
        run%mesh%coordinates(2,node)-run%topology%center(2))
      hole_node(node)=abs(radius-run%config%hole_radius) <= tolerance
    end do
    overlap=count(constrained .and. hole_node)
  end function count_hole_dirichlet_overlap

  pure function metric_values(metrics) result(values)
    type(HoleFieldMetricsType), intent(in) :: metrics
    real(dp) :: values(METRIC_COUNT)
    values=[metrics%total_joule_power,metrics%annulus_joule_mean, &
      metrics%annulus_joule_p95,metrics%annulus_electric_field_mean, &
      metrics%annulus_n_relative_mean,metrics%annulus_temperature_rise_mean]
  end function metric_values

  elemental real(dp) function relative_change(coarse,fine) result(change)
    real(dp), intent(in) :: coarse,fine
    change=abs(fine-coarse)/max(abs(fine),tiny(1.0_dp))
  end function relative_change

  subroutine print_level_metrics(level,run,values)
    integer, intent(in) :: level
    type(G10CaseRunType), intent(in) :: run
    real(dp), intent(in) :: values(METRIC_COUNT)
    character(len=8) :: label
    integer :: metric

    write(label,'("l",i0)') level
    write(*,'(a,i0)') 'g10_'//trim(label)//'_node_count=',run%mesh%number_of_nodes()
    write(*,'(a,i0)') 'g10_'//trim(label)//'_element_count=',run%mesh%number_of_elements()
    write(*,'(a,i0)') 'g10_'//trim(label)//'_hole_faces=',run%topology%hole_face_count
    write(*,'(a,i0)') 'g10_'//trim(label)//'_hole_nodes=',run%topology%hole_node_count
    write(*,'(a,i0)') 'g10_'//trim(label)//'_accepted_points=',run%report%accepted_points
    write(*,'(a,i0)') 'g10_'//trim(label)//'_rejected_attempts=',run%report%rejected_attempts
    call print_metric('g10_'//trim(label)//'_minimum_determinant',run%minimum_determinant)
    call print_metric('g10_'//trim(label)//'_volume_relative_error',run%volume_error)
    call print_metric('g10_'//trim(label)//'_annulus_volume',run%metrics%annulus_volume)
    call print_metric('g10_'//trim(label)//'_maximum_balance', &
      run%report%maximum_balance_defect)
    do metric=1,METRIC_COUNT
      call print_metric('g10_'//trim(label)//'_'//trim(metric_name(metric)),values(metric))
    end do
  end subroutine print_level_metrics

  pure function metric_name(metric) result(name)
    integer, intent(in) :: metric
    character(len=32) :: name
    select case(metric)
    case(1); name='total_joule_power'
    case(2); name='annulus_joule_mean'
    case(3); name='annulus_joule_p95'
    case(4); name='annulus_electric_mean'
    case(5); name='annulus_n_relative_mean'
    case(6); name='annulus_temperature_mean'
    case default; name='unknown'
    end select
  end function metric_name

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

  function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

  subroutine copy_failure(prefix,source,destination)
    character(len=*), intent(in) :: prefix
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,STATUS_VERIFICATION_FAILED,trim(prefix)//trim(source%message))
  end subroutine copy_failure

end module mesh_convergence_check_m

module matched_control_check_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use input_schema_m, only : CaseConfigType,PROBLEM_PLAIN_PLATE, &
    PROBLEM_CENTER_HOLE_PLATE
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType,StateType,LoadStateType
  use status_m, only : SolverStatusType,clear_status,set_status, &
    STATUS_VERIFICATION_FAILED
  use native_input_m, only : read_native_case
  use mesh_checks_m, only : validate_mesh
  use boundary_topology_m, only : BoundaryTopologyAuditType,audit_boundary_topology
  use sparse_pattern_m, only : SparseAssemblyContextType,initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use reference_bar_problem_m, only : initialize_voltage_driven_problem
  use continuation_solver_m, only : ContinuationReportType,solve_load_continuation
  use field_metrics_m, only : HoleFieldMetricsType,evaluate_hole_field_metrics
  implicit none
  private

  integer, parameter :: METRIC_COUNT=9
  integer, parameter :: ANNULUS_SUBDIVISIONS=8
  real(dp), parameter :: VIRTUAL_HOLE_RADIUS=90.0e-9_dp
  real(dp), parameter :: ANNULUS_WIDTH=30.0e-9_dp

  type :: G11CaseRunType
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
  end type G11CaseRunType

  public :: verify_matched_control_g11

contains

  subroutine verify_matched_control_g11(plain_case_path,hole_case_path,status)
    character(len=*), intent(in) :: plain_case_path,hole_case_path
    type(SolverStatusType), intent(out) :: status
    type(G11CaseRunType) :: plain,hole
    type(SolverStatusType) :: local_status
    real(dp) :: plain_values(METRIC_COUNT),hole_values(METRIC_COUNT)
    real(dp) :: contrast(METRIC_COUNT),boundary_mismatch
    real(dp) :: plain_exact_volume,hole_exact_volume
    real(dp) :: exact_annulus_volume,plain_annulus_volume_error,hole_annulus_volume_error
    logical :: passed
    integer :: metric,hole_overlap

    call clear_status(status)
    call prepare_solve_and_measure(plain_case_path,plain,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G11 H00 failed: ',local_status,status); return
    end if
    call prepare_solve_and_measure(hole_case_path,hole,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G11 H02 failed: ',local_status,status); return
    end if
    if (trim(plain%config%problem_type) /= PROBLEM_PLAIN_PLATE .or. &
        trim(hole%config%problem_type) /= PROBLEM_CENTER_HOLE_PLATE) then
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G11 paths must identify H00 plain and H02 center-hole cases.')
      return
    end if

    plain_exact_volume=product(plain%topology%bounds_max-plain%topology%bounds_min)
    hole_exact_volume=((hole%topology%bounds_max(1)-hole%topology%bounds_min(1))* &
      (hole%topology%bounds_max(2)-hole%topology%bounds_min(2))- &
      acos(-1.0_dp)*hole%config%hole_radius**2)* &
      (hole%topology%bounds_max(3)-hole%topology%bounds_min(3))
    plain%volume_error=abs(plain%volume-plain_exact_volume)/plain_exact_volume
    hole%volume_error=abs(hole%volume-hole_exact_volume)/hole_exact_volume
    exact_annulus_volume=acos(-1.0_dp)*((VIRTUAL_HOLE_RADIUS+ANNULUS_WIDTH)**2- &
      VIRTUAL_HOLE_RADIUS**2)*(plain%topology%bounds_max(3)-plain%topology%bounds_min(3))
    plain_annulus_volume_error=abs(plain%metrics%annulus_volume-exact_annulus_volume)/ &
      exact_annulus_volume
    hole_annulus_volume_error=abs(hole%metrics%annulus_volume-exact_annulus_volume)/ &
      exact_annulus_volume
    boundary_mismatch=maximum_outer_boundary_mismatch(plain,hole)
    hole_overlap=count_hole_dirichlet_overlap(hole)
    plain_values=metric_values(plain%metrics)
    hole_values=metric_values(hole%metrics)
    contrast=(hole_values-plain_values)/plain_values

    passed=configurations_match(plain%config,hole%config) .and. &
      plain_topology_passes(plain) .and. hole_topology_passes(hole) .and. &
      plain%minimum_determinant > 0.0_dp .and. hole%minimum_determinant > 0.0_dp .and. &
      plain%volume_error < 1.0e-3_dp .and. hole%volume_error < 1.0e-3_dp .and. &
      plain_annulus_volume_error < 0.02_dp .and. &
      hole_annulus_volume_error < 0.02_dp .and. &
      boundary_mismatch < 1.0e-18_dp .and. hole_overlap == 0 .and. &
      continuation_passes(plain) .and. continuation_passes(hole) .and. &
      all(ieee_is_finite(plain_values)) .and. all(ieee_is_finite(hole_values)) .and. &
      all(plain_values > tiny(1.0_dp)) .and. all(hole_values >= 0.0_dp) .and. &
      all(ieee_is_finite(contrast))

    write(*,'(a)') 'verification_id=G11_matched_control'
    call print_metric('g11_virtual_hole_radius',VIRTUAL_HOLE_RADIUS)
    call print_metric('g11_annulus_width',ANNULUS_WIDTH)
    write(*,'(a,i0)') 'g11_annulus_subdivisions=',ANNULUS_SUBDIVISIONS
    call print_metric('g11_outer_boundary_mismatch',boundary_mismatch)
    call print_metric('g11_exact_annulus_volume',exact_annulus_volume)
    call print_metric('g11_h00_annulus_volume_relative_error',plain_annulus_volume_error)
    call print_metric('g11_h02_annulus_volume_relative_error',hole_annulus_volume_error)
    call print_case('h00',plain)
    call print_case('h02',hole)
    write(*,'(a,i0)') 'g11_h02_hole_dirichlet_overlap=',hole_overlap
    do metric=1,METRIC_COUNT
      call print_metric('g11_h00_'//trim(metric_name(metric)),plain_values(metric))
      call print_metric('g11_h02_'//trim(metric_name(metric)),hole_values(metric))
      call print_metric('g11_contrast_'//trim(metric_name(metric)),contrast(metric))
    end do
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G11 matching, topology, nonlinear, conservation, or metric validity failed.')
    end if
  end subroutine verify_matched_control_g11

  subroutine prepare_solve_and_measure(case_path,run,status)
    character(len=*), intent(in) :: case_path
    type(G11CaseRunType), intent(out) :: run
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
      VIRTUAL_HOLE_RADIUS,ANNULUS_WIDTH,run%metrics,status,ANNULUS_SUBDIVISIONS)
  end subroutine prepare_solve_and_measure

  logical function configurations_match(a,b) result(matches)
    type(CaseConfigType), intent(in) :: a,b
    matches=trim(a%model_level) == trim(b%model_level) .and. &
      trim(a%material_profile) == trim(b%material_profile) .and. &
      a%target_voltage == b%target_voltage .and. &
      a%target_voltage_factor == b%target_voltage_factor .and. &
      a%target_joule_factor == b%target_joule_factor .and. &
      all(a%field_scales == b%field_scales) .and. &
      a%nonlinear_options%max_newton == b%nonlinear_options%max_newton .and. &
      a%nonlinear_options%residual_tolerance == b%nonlinear_options%residual_tolerance .and. &
      a%nonlinear_options%update_tolerance == b%nonlinear_options%update_tolerance .and. &
      a%nonlinear_options%armijo_c1 == b%nonlinear_options%armijo_c1 .and. &
      a%nonlinear_options%backtrack_factor == b%nonlinear_options%backtrack_factor .and. &
      a%nonlinear_options%minimum_step == b%nonlinear_options%minimum_step .and. &
      a%nonlinear_options%maximum_log_density_update == &
        b%nonlinear_options%maximum_log_density_update .and. &
      a%continuation_options%initial_step == b%continuation_options%initial_step .and. &
      a%continuation_options%minimum_step == b%continuation_options%minimum_step .and. &
      a%continuation_options%maximum_step == b%continuation_options%maximum_step .and. &
      a%continuation_options%cutback_factor == b%continuation_options%cutback_factor .and. &
      a%continuation_options%growth_factor == b%continuation_options%growth_factor .and. &
      a%continuation_options%maximum_attempts == b%continuation_options%maximum_attempts .and. &
      a%continuation_options%maximum_newton_per_attempt == &
        b%continuation_options%maximum_newton_per_attempt .and. &
      a%continuation_options%fast_iteration_threshold == &
        b%continuation_options%fast_iteration_threshold .and. &
      a%continuation_options%slow_iteration_threshold == &
        b%continuation_options%slow_iteration_threshold
  end function configurations_match

  logical function plain_topology_passes(run) result(passed)
    type(G11CaseRunType), intent(in) :: run
    passed=run%mesh%number_of_nodes() == 531 .and. &
      run%mesh%number_of_elements() == 64 .and. &
      run%topology%left_face_count == 8 .and. run%topology%right_face_count == 8 .and. &
      run%topology%front_face_count == 8 .and. run%topology%back_face_count == 8 .and. &
      run%topology%top_face_count == 64 .and. run%topology%bottom_face_count == 64 .and. &
      run%topology%hole_face_count == 0 .and. run%topology%hole_node_count == 0 .and. &
      run%topology%unknown_face_count == 0 .and. &
      run%topology%left_node_count == 43 .and. run%topology%right_node_count == 43
  end function plain_topology_passes

  logical function hole_topology_passes(run) result(passed)
    type(G11CaseRunType), intent(in) :: run
    passed=run%mesh%number_of_nodes() == 1056 .and. &
      run%mesh%number_of_elements() == 128 .and. &
      run%topology%left_face_count == 8 .and. run%topology%right_face_count == 8 .and. &
      run%topology%front_face_count == 8 .and. run%topology%back_face_count == 8 .and. &
      run%topology%top_face_count == 128 .and. run%topology%bottom_face_count == 128 .and. &
      run%topology%hole_face_count == 32 .and. run%topology%hole_node_count == 160 .and. &
      run%topology%unknown_face_count == 0 .and. &
      run%topology%left_node_count == 43 .and. run%topology%right_node_count == 43 .and. &
      run%topology%maximum_hole_radius_error <= run%topology%classification_tolerance
  end function hole_topology_passes

  logical function continuation_passes(run) result(passed)
    type(G11CaseRunType), intent(in) :: run
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

  real(dp) function maximum_outer_boundary_mismatch(a,b) result(mismatch)
    type(G11CaseRunType), intent(in) :: a,b
    real(dp) :: local_mismatch
    integer :: axis
    mismatch=0.0_dp
    do axis=1,2
      call compare_side_sets(a%mesh,a%topology,a%topology%bounds_min(axis), &
        b%mesh,b%topology,b%topology%bounds_min(axis),axis,local_mismatch)
      mismatch=max(mismatch,local_mismatch)
      call compare_side_sets(a%mesh,a%topology,a%topology%bounds_max(axis), &
        b%mesh,b%topology,b%topology%bounds_max(axis),axis,local_mismatch)
      mismatch=max(mismatch,local_mismatch)
    end do
  end function maximum_outer_boundary_mismatch

  subroutine compare_side_sets(mesh_a,topology_a,value_a,mesh_b,topology_b, &
                               value_b,axis,mismatch)
    type(MeshType), intent(in) :: mesh_a,mesh_b
    type(BoundaryTopologyAuditType), intent(in) :: topology_a,topology_b
    real(dp), intent(in) :: value_a,value_b
    integer, intent(in) :: axis
    real(dp), intent(out) :: mismatch
    real(dp), allocatable :: points_a(:,:),points_b(:,:)
    real(dp) :: nearest
    integer :: i,j

    call collect_side_points(mesh_a,axis,value_a,topology_a%classification_tolerance,points_a)
    call collect_side_points(mesh_b,axis,value_b,topology_b%classification_tolerance,points_b)
    if (size(points_a,2) /= size(points_b,2)) then
      mismatch=huge(1.0_dp); return
    end if
    mismatch=0.0_dp
    do i=1,size(points_a,2)
      nearest=huge(1.0_dp)
      do j=1,size(points_b,2)
        nearest=min(nearest,sqrt(sum((points_a(:,i)-points_b(:,j))**2)))
      end do
      mismatch=max(mismatch,nearest)
    end do
    do i=1,size(points_b,2)
      nearest=huge(1.0_dp)
      do j=1,size(points_a,2)
        nearest=min(nearest,sqrt(sum((points_b(:,i)-points_a(:,j))**2)))
      end do
      mismatch=max(mismatch,nearest)
    end do
  end subroutine compare_side_sets

  subroutine collect_side_points(mesh,axis,value,tolerance,points)
    type(MeshType), intent(in) :: mesh
    integer, intent(in) :: axis
    real(dp), intent(in) :: value,tolerance
    real(dp), allocatable, intent(out) :: points(:,:)
    integer :: node,index
    allocate(points(3,count(abs(mesh%coordinates(axis,:)-value) <= tolerance)))
    index=0
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(axis,node)-value) > tolerance) cycle
      index=index+1; points(:,index)=mesh%coordinates(:,node)
    end do
  end subroutine collect_side_points

  integer function count_hole_dirichlet_overlap(run) result(overlap)
    type(G11CaseRunType), intent(in) :: run
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
    values=[metrics%total_joule_power,metrics%total_joule_power/metrics%volume, &
      metrics%annulus_joule_mean,metrics%annulus_joule_p95, &
      metrics%annulus_electric_field_mean,metrics%annulus_n_relative_mean, &
      metrics%annulus_p_relative_mean,metrics%annulus_temperature_rise_mean, &
      metrics%annulus_displacement_mean]
  end function metric_values

  pure function metric_name(metric) result(name)
    integer, intent(in) :: metric
    character(len=32) :: name
    select case(metric)
    case(1); name='total_joule_power'
    case(2); name='global_joule_mean'
    case(3); name='annulus_joule_mean'
    case(4); name='annulus_joule_p95'
    case(5); name='annulus_electric_mean'
    case(6); name='annulus_n_relative_mean'
    case(7); name='annulus_p_relative_mean'
    case(8); name='annulus_temperature_rise_mean'
    case(9); name='annulus_displacement_mean'
    case default; name='unknown'
    end select
  end function metric_name

  subroutine print_case(label,run)
    character(len=*), intent(in) :: label
    type(G11CaseRunType), intent(in) :: run
    write(*,'(a,i0)') 'g11_'//trim(label)//'_node_count=',run%mesh%number_of_nodes()
    write(*,'(a,i0)') 'g11_'//trim(label)//'_element_count=',run%mesh%number_of_elements()
    write(*,'(a,i0)') 'g11_'//trim(label)//'_accepted_points=',run%report%accepted_points
    write(*,'(a,i0)') 'g11_'//trim(label)//'_rejected_attempts=',run%report%rejected_attempts
    call print_metric('g11_'//trim(label)//'_volume_relative_error',run%volume_error)
    call print_metric('g11_'//trim(label)//'_maximum_balance', &
      run%report%maximum_balance_defect)
  end subroutine print_case

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

end module matched_control_check_m

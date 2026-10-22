module radius_matrix_check_m
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

  integer, parameter :: RADIUS_COUNT=3,METRIC_COUNT=9
  integer, parameter :: ANNULUS_SUBDIVISIONS=8
  real(dp), parameter :: RADII(RADIUS_COUNT)=[60.0e-9_dp,90.0e-9_dp,120.0e-9_dp]
  real(dp), parameter :: ANNULUS_WIDTH=30.0e-9_dp

  type :: G12RunType
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
    real(dp) :: annulus_volume_error=huge(1.0_dp)
    real(dp) :: boundary_mismatch=huge(1.0_dp)
    integer :: hole_dirichlet_overlap=0
  end type G12RunType

  public :: verify_radius_matrix_g12,verify_radius_controls_g12c

contains

  subroutine verify_radius_matrix_g12(plain_case_path,hole_case_paths,status)
    character(len=*), intent(in) :: plain_case_path
    character(len=*), intent(in) :: hole_case_paths(RADIUS_COUNT)
    type(SolverStatusType), intent(out) :: status
    type(G12RunType) :: plain,hole(RADIUS_COUNT)
    type(SolverStatusType) :: local_status
    real(dp) :: plain_values(METRIC_COUNT,RADIUS_COUNT),hole_values(METRIC_COUNT,RADIUS_COUNT)
    real(dp) :: plain_annulus_errors(RADIUS_COUNT)
    real(dp) :: contrasts(METRIC_COUNT,RADIUS_COUNT),plain_exact_volume
    real(dp) :: exact_volume,exact_annulus_volume
    logical :: passed
    integer :: radius_index,metric

    call clear_status(status)
    call prepare_solve_and_measure(plain_case_path,0.0_dp,plain,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G12 H00 failed: ',local_status,status); return
    end if
    do radius_index=1,RADIUS_COUNT
      call prepare_solve_and_measure(hole_case_paths(radius_index),RADII(radius_index), &
        hole(radius_index),local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure('G12 R'//radius_label(radius_index)//' failed: ',local_status,status)
        return
      end if
    end do

    plain_exact_volume=product(plain%topology%bounds_max-plain%topology%bounds_min)
    plain%volume_error=abs(plain%volume-plain_exact_volume)/plain_exact_volume
    plain%annulus_volume_error=0.0_dp
    plain%boundary_mismatch=0.0_dp
    do radius_index=1,RADIUS_COUNT
      exact_volume=((hole(radius_index)%topology%bounds_max(1)- &
        hole(radius_index)%topology%bounds_min(1))* &
        (hole(radius_index)%topology%bounds_max(2)-hole(radius_index)%topology%bounds_min(2))- &
        acos(-1.0_dp)*hole(radius_index)%config%hole_radius**2)* &
        (hole(radius_index)%topology%bounds_max(3)-hole(radius_index)%topology%bounds_min(3))
      hole(radius_index)%volume_error=abs(hole(radius_index)%volume-exact_volume)/exact_volume
      exact_annulus_volume=acos(-1.0_dp)*((RADII(radius_index)+ANNULUS_WIDTH)**2- &
        RADII(radius_index)**2)*(plain%topology%bounds_max(3)-plain%topology%bounds_min(3))
      hole(radius_index)%annulus_volume_error=abs(hole(radius_index)%metrics%annulus_volume- &
        exact_annulus_volume)/exact_annulus_volume
      call evaluate_hole_field_metrics(plain%mesh,plain%materials,plain%state,plain%load, &
        RADII(radius_index),ANNULUS_WIDTH,plain%metrics,local_status,ANNULUS_SUBDIVISIONS)
      if (.not. local_status%is_ok()) then
        call copy_failure('G12 H00 virtual annulus failed: ',local_status,status); return
      end if
      exact_annulus_volume=acos(-1.0_dp)*((RADII(radius_index)+ANNULUS_WIDTH)**2- &
        RADII(radius_index)**2)*(plain%topology%bounds_max(3)-plain%topology%bounds_min(3))
      plain_annulus_errors(radius_index)=abs(plain%metrics%annulus_volume-exact_annulus_volume)/ &
        exact_annulus_volume
      plain_values(:,radius_index)=metric_values(plain%metrics)
      hole_values(:,radius_index)=metric_values(hole(radius_index)%metrics)
      hole(radius_index)%boundary_mismatch=maximum_outer_boundary_mismatch(plain,hole(radius_index))
    end do
    ! Each hole radius has its own H00 virtual annulus. Validate every
    ! denominator before division and retain every baseline for reporting.
    if (.not. radius_controls_pass(plain_values,plain_annulus_errors)) then
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G12 H00 virtual-annulus volume or denominator invalid at one or more radii.')
      return
    end if
    contrasts=(hole_values-plain_values)/plain_values
    ! The common H00 summary reports the worst error, not the final radius.
    plain%annulus_volume_error=maxval(plain_annulus_errors)

    passed=plain_topology_passes(plain) .and. continuation_passes(plain) .and. &
      plain%minimum_determinant > 0.0_dp .and. plain%volume_error < 1.0e-3_dp .and. &
      radius_controls_pass(plain_values,plain_annulus_errors)
    do radius_index=1,RADIUS_COUNT
      passed=passed .and. hole_topology_passes(hole(radius_index),radius_index) .and. &
        continuation_passes(hole(radius_index)) .and. hole(radius_index)%minimum_determinant > 0.0_dp .and. &
        hole(radius_index)%volume_error < 1.0e-3_dp .and. &
        hole(radius_index)%annulus_volume_error < 0.02_dp .and. &
        hole(radius_index)%boundary_mismatch < 1.0e-18_dp .and. &
        hole(radius_index)%hole_dirichlet_overlap == 0 .and. &
        all(ieee_is_finite(hole_values(:,radius_index))) .and. &
        all(hole_values(:,radius_index) >= 0.0_dp) .and. all(ieee_is_finite(contrasts(:,radius_index)))
    end do

    write(*,'(a)') 'verification_id=G12_radius_matrix'
    call print_metric('g12_annulus_width',ANNULUS_WIDTH)
    write(*,'(a,i0)') 'g12_annulus_subdivisions=',ANNULUS_SUBDIVISIONS
    call print_case('h00',plain)
    do radius_index=1,RADIUS_COUNT
      call print_case('r'//radius_label(radius_index),hole(radius_index))
      call print_metric('g12_r'//radius_label(radius_index)//'_outer_boundary_mismatch', &
        hole(radius_index)%boundary_mismatch)
      call print_metric('g12_r'//radius_label(radius_index)//'_h00_annulus_volume_relative_error', &
        plain_annulus_errors(radius_index))
      do metric=1,METRIC_COUNT
        call print_metric('g12_r'//radius_label(radius_index)//'_h00_'//trim(metric_name(metric)), &
          plain_values(metric,radius_index))
        call print_metric('g12_r'//radius_label(radius_index)//'_h02_'//trim(metric_name(metric)), &
          hole_values(metric,radius_index))
        call print_metric('g12_r'//radius_label(radius_index)//'_contrast_'//trim(metric_name(metric)), &
          contrasts(metric,radius_index))
      end do
    end do
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G12 radius matrix topology, nonlinear, conservation, annulus, or metric validity failed.')
    end if
  end subroutine verify_radius_matrix_g12

  pure logical function radius_controls_pass(values,annulus_errors) result(passed)
    real(dp), intent(in) :: values(METRIC_COUNT,RADIUS_COUNT),annulus_errors(RADIUS_COUNT)
    passed=all(ieee_is_finite(values)) .and. all(values > tiny(1.0_dp)) .and. &
      all(ieee_is_finite(annulus_errors)) .and. all(annulus_errors >= 0.0_dp) .and. &
      all(annulus_errors < 0.02_dp)
  end function radius_controls_pass

  subroutine verify_radius_controls_g12c(status)
    use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan
    type(SolverStatusType), intent(out) :: status
    real(dp) :: values(METRIC_COUNT,RADIUS_COUNT),errors(RADIUS_COUNT)
    real(dp) :: test_values(METRIC_COUNT,RADIUS_COUNT),test_errors(RADIUS_COUNT)
    integer :: radius_index,metric,negative_checks
    logical :: passed

    call clear_status(status)
    values=1.0_dp;errors=0.005_dp
    passed=radius_controls_pass(values,errors)
    negative_checks=0
    do radius_index=1,RADIUS_COUNT
      do metric=1,METRIC_COUNT
        test_values=values;test_values(metric,radius_index)=0.0_dp
        passed=passed .and. .not. radius_controls_pass(test_values,errors)
        test_values(metric,radius_index)=-1.0_dp
        passed=passed .and. .not. radius_controls_pass(test_values,errors)
        test_values(metric,radius_index)=ieee_value(0.0_dp,ieee_quiet_nan)
        passed=passed .and. .not. radius_controls_pass(test_values,errors)
        negative_checks=negative_checks+3
      end do
      test_errors=errors;test_errors(radius_index)=0.02_dp
      passed=passed .and. .not. radius_controls_pass(values,test_errors)
      test_errors(radius_index)=-1.0_dp
      passed=passed .and. .not. radius_controls_pass(values,test_errors)
      test_errors(radius_index)=ieee_value(0.0_dp,ieee_quiet_nan)
      passed=passed .and. .not. radius_controls_pass(values,test_errors)
      negative_checks=negative_checks+3
    end do
    write(*,'(a)') 'verification_id=G12C_radius_control_guard'
    write(*,'(a,i0)') 'g12c_negative_checks=',negative_checks
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'G12C per-radius control guard regression failed.')
    end if
  end subroutine verify_radius_controls_g12c

  subroutine prepare_solve_and_measure(case_path,expected_radius,run,status)
    character(len=*), intent(in) :: case_path
    real(dp), intent(in) :: expected_radius
    type(G12RunType), intent(out) :: run
    type(SolverStatusType), intent(out) :: status
    type(LoadStateType) :: target_load

    call read_native_case(case_path,run%config,run%mesh,status)
    if (.not. status%is_ok()) return
    if (expected_radius > 0.0_dp .and. abs(run%config%hole_radius-expected_radius) > 1.0e-18_dp) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'G12 case radius does not match frozen matrix.')
      return
    end if
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
    if (expected_radius > 0.0_dp) then
      call evaluate_hole_field_metrics(run%mesh,run%materials,run%state,run%load,expected_radius, &
        ANNULUS_WIDTH,run%metrics,status,ANNULUS_SUBDIVISIONS)
    end if
  end subroutine prepare_solve_and_measure

  logical function plain_topology_passes(run) result(passed)
    type(G12RunType), intent(in) :: run
    passed=trim(run%config%problem_type) == PROBLEM_PLAIN_PLATE .and. &
      run%mesh%number_of_nodes() == 531 .and. run%mesh%number_of_elements() == 64 .and. &
      run%topology%left_face_count == 8 .and. run%topology%right_face_count == 8 .and. &
      run%topology%front_face_count == 8 .and. run%topology%back_face_count == 8 .and. &
      run%topology%top_face_count == 64 .and. run%topology%bottom_face_count == 64 .and. &
      run%topology%hole_face_count == 0 .and. run%topology%hole_node_count == 0 .and. &
      run%topology%unknown_face_count == 0 .and. run%topology%left_node_count == 43 .and. &
      run%topology%right_node_count == 43
  end function plain_topology_passes

  logical function hole_topology_passes(run,radius_index) result(passed)
    type(G12RunType), intent(in) :: run
    integer, intent(in) :: radius_index
    passed=trim(run%config%problem_type) == PROBLEM_CENTER_HOLE_PLATE .and. &
      run%mesh%number_of_nodes() == 1056 .and. run%mesh%number_of_elements() == 128 .and. &
      run%topology%left_face_count == 8 .and. run%topology%right_face_count == 8 .and. &
      run%topology%front_face_count == 8 .and. run%topology%back_face_count == 8 .and. &
      run%topology%top_face_count == 128 .and. run%topology%bottom_face_count == 128 .and. &
      run%topology%hole_face_count == 32 .and. run%topology%hole_node_count == 160 .and. &
      run%topology%unknown_face_count == 0 .and. run%topology%left_node_count == 43 .and. &
      run%topology%right_node_count == 43 .and. run%topology%maximum_hole_radius_error <= &
      run%topology%classification_tolerance .and. abs(run%config%hole_radius-RADII(radius_index)) <= 1.0e-18_dp
  end function hole_topology_passes

  logical function continuation_passes(run) result(passed)
    type(G12RunType), intent(in) :: run
    integer :: point
    passed=run%report%converged .and. run%load%voltage_factor == 1.0_dp .and. &
      run%load%joule_factor == 1.0_dp .and. run%report%maximum_balance_defect < 1.0e-8_dp .and. &
      run%report%maximum_boundary_error < 1.0e-14_dp
    do point=1,run%report%accepted_points
      passed=passed .and. run%report%residual_rms(point) < 1.0e-9_dp .and. &
        run%report%update_rms(point) < 1.0e-9_dp .and. run%report%minimum_temperature(point) > 0.0_dp .and. &
        run%report%minimum_n(point) > 0.0_dp .and. run%report%minimum_p(point) > 0.0_dp .and. &
        ieee_is_finite(run%report%maximum_temperature(point)) .and. &
        ieee_is_finite(run%report%maximum_n(point)) .and. ieee_is_finite(run%report%maximum_p(point))
    end do
  end function continuation_passes

  real(dp) function maximum_outer_boundary_mismatch(a,b) result(mismatch)
    type(G12RunType), intent(in) :: a,b
    real(dp) :: local_mismatch
    integer :: axis
    mismatch=0.0_dp
    do axis=1,2
      call compare_side_sets(a%mesh,a%topology,a%topology%bounds_min(axis),b%mesh,b%topology, &
        b%topology%bounds_min(axis),axis,local_mismatch)
      mismatch=max(mismatch,local_mismatch)
      call compare_side_sets(a%mesh,a%topology,a%topology%bounds_max(axis),b%mesh,b%topology, &
        b%topology%bounds_max(axis),axis,local_mismatch)
      mismatch=max(mismatch,local_mismatch)
    end do
  end function maximum_outer_boundary_mismatch

  subroutine compare_side_sets(mesh_a,topology_a,value_a,mesh_b,topology_b,value_b,axis,mismatch)
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
      index=index+1;points(:,index)=mesh%coordinates(:,node)
    end do
  end subroutine collect_side_points

  pure function metric_values(metrics) result(values)
    type(HoleFieldMetricsType), intent(in) :: metrics
    real(dp) :: values(METRIC_COUNT)
    values=[metrics%total_joule_power,metrics%total_joule_power/metrics%volume, &
      metrics%annulus_joule_mean,metrics%annulus_joule_p95,metrics%annulus_electric_field_mean, &
      metrics%annulus_n_relative_mean,metrics%annulus_p_relative_mean, &
      metrics%annulus_temperature_rise_mean,metrics%annulus_displacement_mean]
  end function metric_values

  pure function metric_name(metric) result(name)
    integer, intent(in) :: metric
    character(len=32) :: name
    select case(metric)
    case(1);name='total_joule_power'
    case(2);name='global_joule_mean'
    case(3);name='annulus_joule_mean'
    case(4);name='annulus_joule_p95'
    case(5);name='annulus_electric_mean'
    case(6);name='annulus_n_relative_mean'
    case(7);name='annulus_p_relative_mean'
    case(8);name='annulus_temperature_rise_mean'
    case(9);name='annulus_displacement_mean'
    case default;name='unknown'
    end select
  end function metric_name

  function radius_label(index) result(label)
    integer, intent(in) :: index
    character(len=3) :: label
    write(label,'(i3.3)') nint(RADII(index)*1.0e9_dp)
  end function radius_label

  subroutine print_case(label,run)
    character(len=*), intent(in) :: label
    type(G12RunType), intent(in) :: run
    write(*,'(a,i0)') 'g12_'//trim(label)//'_node_count=',run%mesh%number_of_nodes()
    write(*,'(a,i0)') 'g12_'//trim(label)//'_element_count=',run%mesh%number_of_elements()
    write(*,'(a,i0)') 'g12_'//trim(label)//'_accepted_points=',run%report%accepted_points
    write(*,'(a,i0)') 'g12_'//trim(label)//'_rejected_attempts=',run%report%rejected_attempts
    call print_metric('g12_'//trim(label)//'_volume_relative_error',run%volume_error)
    call print_metric('g12_'//trim(label)//'_annulus_volume_relative_error',run%annulus_volume_error)
    call print_metric('g12_'//trim(label)//'_maximum_balance',run%report%maximum_balance_defect)
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

end module radius_matrix_check_m

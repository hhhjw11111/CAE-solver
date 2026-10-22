module center_hole_check_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE
  use input_schema_m, only : CaseConfigType, PROBLEM_PLAIN_PLATE, &
    PROBLEM_CENTER_HOLE_PLATE
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, DiagnosticsType
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_VERIFICATION_FAILED
  use native_input_m, only : read_native_case
  use mesh_checks_m, only : validate_mesh
  use boundary_topology_m, only : BoundaryTopologyAuditType, audit_boundary_topology
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context, csr_matvec
  use dirichlet_m, only : DirichletBCType, build_fixed_dof_mask
  use reference_bar_problem_m, only : initialize_voltage_driven_problem
  use continuation_solver_m, only : ContinuationReportType, solve_load_continuation
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, &
    assemble_global_system
  use scaling_m, only : build_variable_scaling, equilibrate_csr_system
  implicit none
  private

  integer, parameter :: STEP_COUNT=5
  real(dp), parameter :: STEP_VALUES(STEP_COUNT)=[1.0e-2_dp,3.0e-3_dp, &
    1.0e-3_dp,3.0e-4_dp,1.0e-4_dp]

  type :: G09CaseRunType
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    type(MaterialType), allocatable :: materials(:)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(SparseAssemblyContextType) :: context
    type(ContinuationReportType) :: report
    type(BoundaryTopologyAuditType) :: topology
    real(dp) :: minimum_determinant=0.0_dp
    real(dp) :: volume=0.0_dp
  end type G09CaseRunType

  public :: verify_center_hole_g09

contains

  subroutine verify_center_hole_g09(plain_case_path,hole_case_path,status)
    character(len=*), intent(in) :: plain_case_path,hole_case_path
    type(SolverStatusType), intent(out) :: status
    type(G09CaseRunType) :: plain,hole
    type(SolverStatusType) :: local_status
    real(dp) :: plain_exact_volume,hole_exact_volume,plain_volume_error,hole_volume_error
    real(dp) :: tangent_error
    integer :: tangent_stable_steps,hole_dirichlet_overlap
    logical :: topology_pass,solver_pass,passed

    call clear_status(status)
    call prepare_and_solve(plain_case_path,plain,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G09 H00 failed: ',local_status,status); return
    end if
    call prepare_and_solve(hole_case_path,hole,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G09 H02 failed: ',local_status,status); return
    end if
    if (trim(plain%config%problem_type) /= PROBLEM_PLAIN_PLATE .or. &
        trim(hole%config%problem_type) /= PROBLEM_CENTER_HOLE_PLATE) then
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G09 case paths do not identify the frozen H00/H02 problem types.')
      return
    end if

    plain_exact_volume=product(plain%topology%bounds_max-plain%topology%bounds_min)
    hole_exact_volume=((hole%topology%bounds_max(1)-hole%topology%bounds_min(1))* &
      (hole%topology%bounds_max(2)-hole%topology%bounds_min(2))- &
      acos(-1.0_dp)*hole%config%hole_radius**2)* &
      (hole%topology%bounds_max(3)-hole%topology%bounds_min(3))
    plain_volume_error=abs(plain%volume-plain_exact_volume)/plain_exact_volume
    hole_volume_error=abs(hole%volume-hole_exact_volume)/hole_exact_volume
    hole_dirichlet_overlap=count_hole_dirichlet_overlap(hole)
    call audit_final_jacobian(hole,tangent_error,tangent_stable_steps,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure('G09 hole Jacobian audit failed: ',local_status,status); return
    end if

    topology_pass=plain%mesh%number_of_nodes() == 155 .and. &
      plain%mesh%number_of_elements() == 16 .and. &
      plain%topology%left_face_count == 4 .and. plain%topology%right_face_count == 4 .and. &
      plain%topology%front_face_count == 4 .and. plain%topology%back_face_count == 4 .and. &
      plain%topology%bottom_face_count == 16 .and. plain%topology%top_face_count == 16 .and. &
      plain%topology%hole_face_count == 0 .and. plain%topology%hole_node_count == 0 .and. &
      plain%topology%unknown_face_count == 0 .and. &
      plain%topology%left_node_count == 23 .and. plain%topology%right_node_count == 23 .and. &
      hole%mesh%number_of_nodes() == 304 .and. hole%mesh%number_of_elements() == 32 .and. &
      hole%topology%left_face_count == 4 .and. hole%topology%right_face_count == 4 .and. &
      hole%topology%front_face_count == 4 .and. hole%topology%back_face_count == 4 .and. &
      hole%topology%bottom_face_count == 32 .and. hole%topology%top_face_count == 32 .and. &
      hole%topology%hole_face_count == 16 .and. hole%topology%hole_node_count == 80 .and. &
      hole%topology%unknown_face_count == 0 .and. &
      hole%topology%left_node_count == 23 .and. hole%topology%right_node_count == 23 .and. &
      hole%topology%maximum_hole_radius_error <= hole%topology%classification_tolerance .and. &
      hole_dirichlet_overlap == 0
    solver_pass=continuation_passes(plain) .and. continuation_passes(hole)
    passed=topology_pass .and. plain%minimum_determinant > 0.0_dp .and. &
      hole%minimum_determinant > 0.0_dp .and. plain_volume_error < 1.0e-12_dp .and. &
      hole_volume_error < 1.0e-3_dp .and. tangent_error < 1.0e-6_dp .and. &
      tangent_stable_steps >= 3 .and. solver_pass

    write(*,'(a)') 'verification_id=G09_center_hole_geometry'
    call print_case_metrics('h00',plain)
    call print_case_metrics('h02',hole)
    call print_metric('g09_h00_volume_relative_error',plain_volume_error)
    call print_metric('g09_h02_volume_relative_error',hole_volume_error)
    call print_metric('g09_h02_hole_radius_error',hole%topology%maximum_hole_radius_error)
    write(*,'(a,i0)') 'g09_h02_hole_dirichlet_overlap=',hole_dirichlet_overlap
    call print_metric('g09_h02_global_jacobian_error',tangent_error)
    write(*,'(a,i0)') 'g09_h02_global_jacobian_stable_steps=',tangent_stable_steps
    call print_metric('g09_h00_final_max_temperature', &
      plain%report%maximum_temperature(plain%report%accepted_points))
    call print_metric('g09_h02_final_max_temperature', &
      hole%report%maximum_temperature(hole%report%accepted_points))
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G09 geometry, topology, Jacobian, continuation, or conservation requirement failed.')
    end if
  end subroutine verify_center_hole_g09

  subroutine prepare_and_solve(case_path,run,status)
    character(len=*), intent(in) :: case_path
    type(G09CaseRunType), intent(out) :: run
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
  end subroutine prepare_and_solve

  subroutine audit_final_jacobian(run,minimum_error,stable_steps,status)
    type(G09CaseRunType), intent(in) :: run
    real(dp), intent(out) :: minimum_error
    integer, intent(out) :: stable_steps
    type(SolverStatusType), intent(out) :: status
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: base_system,plus_system,minus_system
    type(DiagnosticsType) :: diagnostics
    type(StateType) :: plus_state,minus_state
    logical, allocatable :: fixed(:)
    real(dp), allocatable :: variable_scales(:),residual_scales(:),scaled_values(:)
    real(dp), allocatable :: scaled_residual(:),direction(:),analytic(:),finite_difference(:)
    real(dp) :: direction_norm,error
    integer :: equation,field,step_index

    call clear_status(status)
    request%need_jacobian=.true.
    call assemble_global_system(run%mesh,run%materials,run%state,run%load,run%context, &
      run%boundary_conditions,request,base_system,diagnostics,status)
    if (.not. status%is_ok()) return
    call build_fixed_dof_mask(run%boundary_conditions,run%mesh%number_of_nodes(),fixed,status)
    if (.not. status%is_ok()) return
    call build_variable_scaling(run%mesh%number_of_nodes(),run%config%field_scales, &
      variable_scales,status)
    if (.not. status%is_ok()) return
    call equilibrate_csr_system(run%context,base_system%jacobian_values, &
      base_system%residual,variable_scales,scaled_values,scaled_residual, &
      residual_scales,status)
    if (.not. status%is_ok()) return
    allocate(direction(run%context%number_of_equations), &
      analytic(run%context%number_of_equations), &
      finite_difference(run%context%number_of_equations), &
      plus_state%x(run%context%number_of_equations), &
      minus_state%x(run%context%number_of_equations))
    do equation=1,size(direction)
      field=mod(equation-1,NDOF_NODE)+1
      direction(equation)=run%config%field_scales(field)* &
        (sin(0.13_dp*real(equation,dp))+cos(0.29_dp*real(equation,dp)))
    end do
    where(fixed) direction=0.0_dp
    direction_norm=sqrt(sum((direction/variable_scales)**2))
    direction=direction/direction_norm
    call csr_matvec(run%context,base_system%jacobian_values,direction,analytic,status)
    if (.not. status%is_ok()) return

    minimum_error=huge(1.0_dp); stable_steps=0; request%need_jacobian=.false.
    do step_index=1,STEP_COUNT
      plus_state%x=run%state%x+STEP_VALUES(step_index)*direction
      minus_state%x=run%state%x-STEP_VALUES(step_index)*direction
      call assemble_global_system(run%mesh,run%materials,plus_state,run%load,run%context, &
        run%boundary_conditions,request,plus_system,diagnostics,status)
      if (.not. status%is_ok()) return
      call assemble_global_system(run%mesh,run%materials,minus_state,run%load,run%context, &
        run%boundary_conditions,request,minus_system,diagnostics,status)
      if (.not. status%is_ok()) return
      finite_difference=(plus_system%residual-minus_system%residual)/ &
        (2.0_dp*STEP_VALUES(step_index))
      error=sqrt(sum(((finite_difference-analytic)/residual_scales)**2))/ &
        max(1.0_dp,sqrt(sum((analytic/residual_scales)**2)))
      minimum_error=min(minimum_error,error)
      if (error < 1.0e-5_dp) stable_steps=stable_steps+1
    end do
  end subroutine audit_final_jacobian

  integer function count_hole_dirichlet_overlap(run) result(overlap)
    type(G09CaseRunType), intent(in) :: run
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

  logical function continuation_passes(run) result(passed)
    type(G09CaseRunType), intent(in) :: run
    integer :: point

    passed=run%report%converged .and. run%load%voltage_factor == 1.0_dp .and. &
      run%load%joule_factor == 1.0_dp .and. run%report%maximum_balance_defect < 1.0e-8_dp .and. &
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

  subroutine print_case_metrics(label,run)
    character(len=*), intent(in) :: label
    type(G09CaseRunType), intent(in) :: run
    write(*,'(a,i0)') 'g09_'//trim(label)//'_node_count=',run%mesh%number_of_nodes()
    write(*,'(a,i0)') 'g09_'//trim(label)//'_element_count=',run%mesh%number_of_elements()
    write(*,'(a,i0)') 'g09_'//trim(label)//'_boundary_faces=', &
      run%topology%boundary_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_left_faces=',run%topology%left_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_right_faces=',run%topology%right_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_top_faces=',run%topology%top_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_bottom_faces=',run%topology%bottom_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_hole_faces=',run%topology%hole_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_hole_nodes=',run%topology%hole_node_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_unknown_faces=',run%topology%unknown_face_count
    write(*,'(a,i0)') 'g09_'//trim(label)//'_accepted_points=',run%report%accepted_points
    write(*,'(a,i0)') 'g09_'//trim(label)//'_rejected_attempts=',run%report%rejected_attempts
    call print_metric('g09_'//trim(label)//'_minimum_determinant',run%minimum_determinant)
    call print_metric('g09_'//trim(label)//'_volume',run%volume)
    call print_metric('g09_'//trim(label)//'_maximum_balance', &
      run%report%maximum_balance_defect)
  end subroutine print_case_metrics

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

end module center_hole_check_m

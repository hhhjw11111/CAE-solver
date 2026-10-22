module global_block_jacobian_check_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : NDOF_NODE, F_U1, F_U2, F_U3, F_PHI, F_ETA_N, F_ETA_P, &
    F_THETA, global_dof
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, DiagnosticsType
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context, csr_matvec
  use dirichlet_m, only : DirichletBCType, project_dirichlet_state, build_fixed_dof_mask
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, &
    assemble_global_system
  use scaling_m, only : build_variable_scaling, equilibrate_csr_system
  use material_tangent_check_m, only : define_test_material
  use global_assembly_check_m, only : build_two_element_mesh
  implicit none
  private

  integer, parameter :: ROW_GROUPS=5
  integer, parameter :: INPUT_FIELDS=7
  integer, parameter :: STEP_COUNT=9
  real(dp), parameter :: STEP_VALUES(STEP_COUNT)=[1.0e-2_dp,3.0e-3_dp,1.0e-3_dp, &
    3.0e-4_dp,1.0e-4_dp,3.0e-5_dp,1.0e-5_dp,3.0e-6_dp,1.0e-6_dp]
  real(dp), parameter :: FIELD_SCALES(INPUT_FIELDS)=[1.0e-10_dp,1.0e-10_dp, &
    1.0e-10_dp,1.0e-1_dp,5.0e-2_dp,5.0e-2_dp,5.0_dp]
  character(len=8), parameter :: ROW_NAMES(ROW_GROUPS)=[character(len=8) :: &
    'u','phi','eta_n','eta_p','theta']
  character(len=8), parameter :: FIELD_NAMES(INPUT_FIELDS)=[character(len=8) :: &
    'u1','u2','u3','phi','eta_n','eta_p','theta']

  public :: verify_global_block_jacobian_g06

contains

  subroutine verify_global_block_jacobian_g06(status,supg,directional)
    type(SolverStatusType), intent(out) :: status
    logical, intent(in), optional :: supg
    logical, intent(in), optional :: directional
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: base_state,plus_state,minus_state
    type(LoadStateType) :: load
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: base_system,plus_system,minus_system
    type(DiagnosticsType) :: diagnostics
    type(SolverStatusType) :: local_status
    logical, allocatable :: fixed(:)
    integer, allocatable :: row_indices(:,:)
    integer :: row_sizes(ROW_GROUPS)
    real(dp), allocatable :: direction(:),analytic_action(:),finite_difference(:)
    real(dp), allocatable :: variable_scales(:),residual_scales(:)
    real(dp), allocatable :: scaled_values(:),scaled_residual(:)
    real(dp) :: minimum_error(ROW_GROUPS,INPUT_FIELDS)
    real(dp) :: best_step(ROW_GROUPS,INPUT_FIELDS),physical_error(ROW_GROUPS,ROW_GROUPS)
    real(dp) :: current_error,zero_leakage,maximum_zero_leakage,fixed_direction_max
    integer :: stable_count(ROW_GROUPS,INPUT_FIELDS)
    integer :: row_group,input_field,step_index,node_count,minimum_stable_count
    logical :: passed

    call clear_status(status)
    call build_two_element_mesh(mesh)
    mesh%coordinates=1.0e-6_dp*mesh%coordinates
    call define_test_material(materials(1))
    if(present(supg)) materials(1)%transport_supg=supg
    if(present(directional)) materials(1)%supg_directional_tau=directional
    materials(1)%reference_p=2.0e20_dp
    node_count=mesh%number_of_nodes()
    allocate(base_state%x(NDOF_NODE*node_count))
    call define_micro_state(mesh,base_state)
    load%voltage_factor=0.60_dp
    load%joule_factor=0.77_dp
    call define_micro_boundary_conditions(mesh,boundary_conditions,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call project_dirichlet_state(base_state,boundary_conditions,load,node_count,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    request%need_jacobian=.true.
    call assemble_global_system(mesh,materials,base_state,load,context,boundary_conditions, &
      request,base_system,diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call build_fixed_dof_mask(boundary_conditions,node_count,fixed,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call build_row_groups(node_count,row_indices,row_sizes)
    call build_variable_scaling(node_count,FIELD_SCALES,variable_scales,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call equilibrate_csr_system(context,base_system%jacobian_values,base_system%residual, &
      variable_scales,scaled_values,scaled_residual,residual_scales,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if

    allocate(direction(context%number_of_equations),analytic_action(context%number_of_equations), &
      finite_difference(context%number_of_equations),plus_state%x(context%number_of_equations), &
      minus_state%x(context%number_of_equations))
    minimum_error=huge(1.0_dp)
    best_step=0.0_dp
    stable_count=0
    maximum_zero_leakage=0.0_dp
    fixed_direction_max=0.0_dp
    request%need_jacobian=.false.

    do input_field=1,INPUT_FIELDS
      call build_single_field_direction(input_field,fixed,direction)
      if (any(fixed)) fixed_direction_max=max(fixed_direction_max,maxval(abs(direction),mask=fixed))
      call csr_matvec(context,base_system%jacobian_values,direction,analytic_action,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      do step_index=1,STEP_COUNT
        plus_state%x=base_state%x+STEP_VALUES(step_index)*direction
        minus_state%x=base_state%x-STEP_VALUES(step_index)*direction
        call assemble_global_system(mesh,materials,plus_state,load,context,boundary_conditions, &
          request,plus_system,diagnostics,local_status)
        if (.not. local_status%is_ok()) then
          call copy_failure(local_status,status); return
        end if
        call assemble_global_system(mesh,materials,minus_state,load,context,boundary_conditions, &
          request,minus_system,diagnostics,local_status)
        if (.not. local_status%is_ok()) then
          call copy_failure(local_status,status); return
        end if
        finite_difference=(plus_system%residual-minus_system%residual)/ &
          (2.0_dp*STEP_VALUES(step_index))
        do row_group=1,ROW_GROUPS
          current_error=scaled_block_error(finite_difference,analytic_action,residual_scales, &
            row_indices(1:row_sizes(row_group),row_group))
          if (current_error < minimum_error(row_group,input_field)) then
            minimum_error(row_group,input_field)=current_error
            best_step(row_group,input_field)=STEP_VALUES(step_index)
          end if
          if (current_error < 1.0e-5_dp) then
            stable_count(row_group,input_field)=stable_count(row_group,input_field)+1
          end if
          if (structural_zero(row_group,input_field)) then
            zero_leakage=scaled_action_norm(finite_difference,residual_scales, &
              row_indices(1:row_sizes(row_group),row_group))
            maximum_zero_leakage=max(maximum_zero_leakage,zero_leakage)
          end if
        end do
      end do
    end do

    call collapse_physical_columns(minimum_error,physical_error)
    minimum_stable_count=STEP_COUNT
    do input_field=1,INPUT_FIELDS
      do row_group=1,ROW_GROUPS
        if (.not. structural_zero(row_group,input_field)) then
          minimum_stable_count=min(minimum_stable_count,stable_count(row_group,input_field))
        end if
      end do
    end do
    passed=maxval(minimum_error) < 1.0e-6_dp .and. minimum_stable_count >= 3 .and. &
      maximum_zero_leakage < 1.0e-12_dp .and. fixed_direction_max == 0.0_dp

    write(*,'(a)') 'verification_id=G06_global_block_jacobian'
    write(*,'(a,i0)') 'g06_node_count=',node_count
    write(*,'(a,i0)') 'g06_equation_count=',context%number_of_equations
    do input_field=1,INPUT_FIELDS
      do row_group=1,ROW_GROUPS
        call print_named_metric('g06_error',ROW_NAMES(row_group),FIELD_NAMES(input_field), &
          minimum_error(row_group,input_field))
        call print_named_metric('g06_best_h',ROW_NAMES(row_group),FIELD_NAMES(input_field), &
          best_step(row_group,input_field))
        write(*,'(a,i0)') 'g06_stable_steps_'//trim(ROW_NAMES(row_group))//'_'// &
          trim(FIELD_NAMES(input_field))//'=',stable_count(row_group,input_field)
      end do
    end do
    do row_group=1,ROW_GROUPS
      write(*,'(a,5(1x,es16.8e3))') 'g06_physical_row_'//trim(ROW_NAMES(row_group))//'=', &
        physical_error(row_group,1),physical_error(row_group,2), &
        physical_error(row_group,3),physical_error(row_group,4),physical_error(row_group,5)
    end do
    call print_metric('g06_maximum_block_error',maxval(minimum_error))
    call print_metric('g06_maximum_structural_zero_leakage',maximum_zero_leakage)
    call print_metric('g06_fixed_direction_max',fixed_direction_max)
    write(*,'(a,i0)') 'g06_minimum_active_stable_steps=',minimum_stable_count
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G06 block Jacobian error, stable-step count, or structural-zero leakage failed.')
    end if
  end subroutine verify_global_block_jacobian_g06

  subroutine define_micro_state(mesh,state)
    type(MeshType), intent(in) :: mesh
    type(StateType), intent(inout) :: state
    real(dp) :: x,y,z,xn,yn,zn
    integer :: node

    state%x=0.0_dp
    do node=1,mesh%number_of_nodes()
      x=mesh%coordinates(1,node); y=mesh%coordinates(2,node); z=mesh%coordinates(3,node)
      xn=x/1.0e-6_dp; yn=y/1.0e-6_dp; zn=z/1.0e-6_dp
      state%x(global_dof(node,F_U1))=1.0e-4_dp*x+0.2e-4_dp*y-0.1e-4_dp*z
      state%x(global_dof(node,F_U2))=-0.3e-4_dp*x+0.8e-4_dp*y+0.2e-4_dp*z
      state%x(global_dof(node,F_U3))=0.1e-4_dp*x-0.2e-4_dp*y+0.6e-4_dp*z
      state%x(global_dof(node,F_PHI))=-1.2e5_dp*x+0.8e5_dp*y-0.6e5_dp*z
      state%x(global_dof(node,F_ETA_N))=0.12_dp+0.05_dp*xn-0.03_dp*yn+0.02_dp*zn
      state%x(global_dof(node,F_ETA_P))=-0.18_dp-0.04_dp*xn+0.06_dp*yn-0.01_dp*zn
      state%x(global_dof(node,F_THETA))=45.0_dp+3.0_dp*xn-2.0_dp*yn+1.5_dp*zn
    end do
  end subroutine define_micro_state

  subroutine define_micro_boundary_conditions(mesh,boundary_conditions,status)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), allocatable, intent(out) :: boundary_conditions(:)
    type(SolverStatusType), intent(out) :: status
    integer :: left_node,right_node,boundary

    call clear_status(status)
    left_node=find_node(mesh,0.0_dp,0.0_dp,0.0_dp)
    right_node=find_node(mesh,2.0e-6_dp,0.0_dp,0.0_dp)
    if (left_node == 0 .or. right_node == 0) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'G06 reference boundary node was not found.')
      return
    end if
    allocate(boundary_conditions(8))
    do boundary=1,8
      allocate(boundary_conditions(boundary)%node_ids(1))
    end do
    boundary_conditions(1)%field_id=F_U1; boundary_conditions(1)%node_ids=left_node
    boundary_conditions(2)%field_id=F_U2; boundary_conditions(2)%node_ids=left_node
    boundary_conditions(3)%field_id=F_U3; boundary_conditions(3)%node_ids=left_node
    boundary_conditions(4)%field_id=F_PHI; boundary_conditions(4)%node_ids=left_node
    boundary_conditions(5)%field_id=F_ETA_N; boundary_conditions(5)%node_ids=left_node
    boundary_conditions(5)%value_start=0.12_dp; boundary_conditions(5)%value_end=0.12_dp
    boundary_conditions(6)%field_id=F_ETA_P; boundary_conditions(6)%node_ids=left_node
    boundary_conditions(6)%value_start=-0.18_dp; boundary_conditions(6)%value_end=-0.18_dp
    boundary_conditions(7)%field_id=F_THETA; boundary_conditions(7)%node_ids=left_node
    boundary_conditions(7)%value_start=45.0_dp; boundary_conditions(7)%value_end=45.0_dp
    boundary_conditions(8)%field_id=F_PHI; boundary_conditions(8)%node_ids=right_node
    boundary_conditions(8)%value_start=0.0_dp; boundary_conditions(8)%value_end=-0.24_dp
  end subroutine define_micro_boundary_conditions

  integer function find_node(mesh,target_x,target_y,target_z) result(node_id)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: target_x,target_y,target_z
    integer :: node
    node_id=0
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(1,node)-target_x) <= 1.0e-16_dp .and. &
          abs(mesh%coordinates(2,node)-target_y) <= 1.0e-16_dp .and. &
          abs(mesh%coordinates(3,node)-target_z) <= 1.0e-16_dp) then
        node_id=node
        return
      end if
    end do
  end function find_node

  subroutine build_row_groups(node_count,indices,sizes)
    integer, intent(in) :: node_count
    integer, allocatable, intent(out) :: indices(:,:)
    integer, intent(out) :: sizes(ROW_GROUPS)
    integer :: node,counter,component

    sizes=[3*node_count,node_count,node_count,node_count,node_count]
    allocate(indices(maxval(sizes),ROW_GROUPS))
    indices=0
    counter=0
    do node=1,node_count
      do component=F_U1,F_U3
        counter=counter+1
        indices(counter,1)=global_dof(node,component)
      end do
      indices(node,2)=global_dof(node,F_PHI)
      indices(node,3)=global_dof(node,F_ETA_N)
      indices(node,4)=global_dof(node,F_ETA_P)
      indices(node,5)=global_dof(node,F_THETA)
    end do
  end subroutine build_row_groups

  subroutine build_single_field_direction(field,fixed,direction)
    integer, intent(in) :: field
    logical, intent(in) :: fixed(:)
    real(dp), intent(out) :: direction(:)
    real(dp) :: norm_value
    integer :: equation

    direction=0.0_dp
    do equation=1,size(direction)
      if (mod(equation-1,NDOF_NODE)+1 /= field .or. fixed(equation)) cycle
      direction(equation)=sin(0.17_dp*real(equation,dp))+cos(0.31_dp*real(equation,dp))
    end do
    norm_value=sqrt(sum(direction*direction))
    if (norm_value > 0.0_dp) direction=FIELD_SCALES(field)*direction/norm_value
  end subroutine build_single_field_direction

  real(dp) function scaled_block_error(finite_difference,analytic,residual_scales,indices) &
    result(error)
    real(dp), intent(in) :: finite_difference(:),analytic(:),residual_scales(:)
    integer, intent(in) :: indices(:)
    real(dp) :: difference_norm,analytic_norm
    difference_norm=sqrt(sum(((finite_difference(indices)-analytic(indices))/ &
      residual_scales(indices))**2))
    analytic_norm=sqrt(sum((analytic(indices)/residual_scales(indices))**2))
    error=difference_norm/max(1.0_dp,analytic_norm)
  end function scaled_block_error

  real(dp) function scaled_action_norm(action,residual_scales,indices) result(value)
    real(dp), intent(in) :: action(:),residual_scales(:)
    integer, intent(in) :: indices(:)
    value=sqrt(sum((action(indices)/residual_scales(indices))**2))
  end function scaled_action_norm

  logical function structural_zero(row_group,input_field)
    integer, intent(in) :: row_group,input_field
    select case(row_group)
    case(1)
      structural_zero=input_field == F_ETA_N .or. input_field == F_ETA_P
    case(2)
      structural_zero=.false.
    case(3)
      structural_zero=input_field <= F_U3 .or. input_field == F_ETA_P
    case(4)
      structural_zero=input_field <= F_U3 .or. input_field == F_ETA_N
    case(5)
      structural_zero=input_field <= F_U3
    case default
      structural_zero=.false.
    end select
  end function structural_zero

  subroutine collapse_physical_columns(field_error,physical_error)
    real(dp), intent(in) :: field_error(ROW_GROUPS,INPUT_FIELDS)
    real(dp), intent(out) :: physical_error(ROW_GROUPS,ROW_GROUPS)
    physical_error(:,1)=maxval(field_error(:,F_U1:F_U3),dim=2)
    physical_error(:,2)=field_error(:,F_PHI)
    physical_error(:,3)=field_error(:,F_ETA_N)
    physical_error(:,4)=field_error(:,F_ETA_P)
    physical_error(:,5)=field_error(:,F_THETA)
  end subroutine collapse_physical_columns

  subroutine print_named_metric(prefix,row_name,column_name,value)
    character(len=*), intent(in) :: prefix,row_name,column_name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(prefix)//'_'//trim(row_name)//'_'//trim(column_name)//'=',value
  end subroutine print_named_metric

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

  subroutine copy_failure(source,destination)
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,source%code,trim(source%message))
  end subroutine copy_failure

end module global_block_jacobian_check_m

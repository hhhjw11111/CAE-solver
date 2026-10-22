module global_assembly_check_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : NDOF_NODE, F_U1, F_U2, F_U3, F_PHI, F_ETA_N, F_ETA_P, &
    F_THETA, global_dof, local_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use hex20_shape_m, only : hex20_natural_coordinates
  use dof_map_m, only : HEX20_NDOF
  use model_types_m, only : MaterialType, StateType, LoadStateType, DiagnosticsType
  use element_kernel_m, only : ElementRequestType, ElementDiagnosticsType, evaluate_hex20_element
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context, &
    csr_matvec, find_csr_position
  use dirichlet_m, only : DirichletBCType, project_dirichlet_state, build_fixed_dof_mask
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, &
    assemble_global_system
  use material_tangent_check_m, only : define_test_material
  implicit none
  private

  integer, parameter :: NUMBER_OF_GROUPS = 5
  public :: verify_global_g04, build_two_element_mesh, define_global_state

contains

  subroutine verify_global_g04(status)
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: base_state, plus_state, minus_state
    type(LoadStateType) :: load
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: base_system, plus_system, minus_system, unconstrained_system
    type(DiagnosticsType) :: diagnostics
    type(SolverStatusType) :: local_status
    logical, allocatable :: fixed(:)
    real(dp), allocatable :: direction(:), raw(:), analytic_action(:), finite_difference(:)
    real(dp) :: minimum_error(NUMBER_OF_GROUPS), current_error(NUMBER_OF_GROUPS)
    real(dp) :: step_values(8), step, csr_structure_error, dirichlet_row_error
    real(dp) :: shared_scatter_error, electron_balance, hole_balance
    integer, allocatable :: indices_u(:), indices_phi(:), indices_n(:), indices_p(:), indices_T(:)
    integer :: equation, field, step_index, node_count
    logical :: passed

    call clear_status(status)
    call build_two_element_mesh(mesh)
    call define_test_material(materials(1))
    materials(1)%reference_p = 2.0e20_dp
    node_count = mesh%number_of_nodes()
    allocate(base_state%x(NDOF_NODE*node_count))
    call define_global_state(mesh,base_state)
    load%voltage_factor = 0.60_dp
    load%joule_factor = 0.77_dp
    call define_boundary_conditions(mesh,boundary_conditions)
    call project_dirichlet_state(base_state,boundary_conditions,load,node_count,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if

    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    request%need_jacobian = .true.
    call assemble_global_system(mesh,materials,base_state,load,context,boundary_conditions, &
      request,base_system,diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call assemble_global_system(mesh,materials,base_state,load,context,request=request, &
      system=unconstrained_system,diagnostics=diagnostics,status=local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if

    call build_fixed_dof_mask(boundary_conditions,node_count,fixed,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call build_global_field_indices(node_count,indices_u,indices_phi,indices_n,indices_p,indices_T)
    allocate(direction(context%number_of_equations),raw(context%number_of_equations), &
      analytic_action(context%number_of_equations),finite_difference(context%number_of_equations))
    do equation = 1, context%number_of_equations
      raw(equation) = sin(0.17_dp*real(equation,dp)) + cos(0.31_dp*real(equation,dp))
    end do
    raw = raw/sqrt(sum(raw*raw))
    do equation = 1, context%number_of_equations
      field = mod(equation-1,NDOF_NODE)+1
      select case(field)
      case(F_U1,F_U2,F_U3)
        direction(equation) = raw(equation)*1.0e-5_dp
      case(F_PHI)
        direction(equation) = raw(equation)*1.0e3_dp
      case(F_ETA_N,F_ETA_P)
        direction(equation) = raw(equation)*1.0e-2_dp
      case(F_THETA)
        direction(equation) = raw(equation)
      end select
      if (fixed(equation)) direction(equation) = 0.0_dp
    end do
    call csr_matvec(context,base_system%jacobian_values,direction,analytic_action,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if

    allocate(plus_state%x(context%number_of_equations),minus_state%x(context%number_of_equations))
    step_values = [1.0e-2_dp,3.0e-3_dp,1.0e-3_dp,3.0e-4_dp,1.0e-4_dp, &
                   3.0e-5_dp,1.0e-5_dp,3.0e-6_dp]
    minimum_error = huge(1.0_dp)
    request%need_jacobian = .false.
    do step_index = 1, size(step_values)
      step = step_values(step_index)
      plus_state%x = base_state%x+step*direction
      minus_state%x = base_state%x-step*direction
      call assemble_global_system(mesh,materials,plus_state,load,context,boundary_conditions, &
        request,plus_system,diagnostics,local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
        return
      end if
      call assemble_global_system(mesh,materials,minus_state,load,context,boundary_conditions, &
        request,minus_system,diagnostics,local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
        return
      end if
      finite_difference = (plus_system%residual-minus_system%residual)/(2.0_dp*step)
      current_error(1) = relative_error(finite_difference(indices_u),analytic_action(indices_u))
      current_error(2) = relative_error(finite_difference(indices_phi),analytic_action(indices_phi))
      current_error(3) = relative_error(finite_difference(indices_n),analytic_action(indices_n))
      current_error(4) = relative_error(finite_difference(indices_p),analytic_action(indices_p))
      current_error(5) = relative_error(finite_difference(indices_T),analytic_action(indices_T))
      minimum_error = min(minimum_error,current_error)
    end do

    csr_structure_error = check_csr_structure(context)
    dirichlet_row_error = check_dirichlet_rows(context,base_system,fixed)
    shared_scatter_error = check_shared_node_scatter(mesh,materials(1),base_state,load, &
      context,unconstrained_system,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    electron_balance = abs(sum(unconstrained_system%residual(indices_n))) / &
      max(1.0_dp,vector_norm(unconstrained_system%residual(indices_n)))
    hole_balance = abs(sum(unconstrained_system%residual(indices_p))) / &
      max(1.0_dp,vector_norm(unconstrained_system%residual(indices_p)))

    passed = node_count == 32 .and. context%number_of_equations == 224 .and. &
      maxval(minimum_error) < 8.0e-7_dp .and. csr_structure_error == 0.0_dp .and. &
      dirichlet_row_error < 1.0e-14_dp .and. shared_scatter_error < 1.0e-13_dp .and. &
      electron_balance < 5.0e-12_dp .and. hole_balance < 5.0e-12_dp

    write(*,'(a)') 'verification_id=G04_global_assembly'
    write(*,'(a,i0)') 'global_node_count=',node_count
    write(*,'(a,i0)') 'global_equation_count=',context%number_of_equations
    write(*,'(a,i0)') 'global_csr_nnz=',context%number_of_nonzeros
    call print_metric('mechanical_row_tangent_error',minimum_error(1))
    call print_metric('potential_row_tangent_error',minimum_error(2))
    call print_metric('electron_row_tangent_error',minimum_error(3))
    call print_metric('hole_row_tangent_error',minimum_error(4))
    call print_metric('thermal_row_tangent_error',minimum_error(5))
    call print_metric('csr_structure_error',csr_structure_error)
    call print_metric('dirichlet_row_error',dirichlet_row_error)
    call print_metric('shared_node_scatter_error',shared_scatter_error)
    call print_metric('global_electron_balance',electron_balance)
    call print_metric('global_hole_balance',hole_balance)
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G04 CSR, scatter, Dirichlet, global tangent, or conservation check failed.')
    end if
  end subroutine verify_global_g04

  subroutine build_two_element_mesh(mesh)
    type(MeshType), intent(out) :: mesh
    real(dp) :: natural_nodes(3,HEX20_NNODE), temporary_coordinates(3,2*HEX20_NNODE)
    real(dp) :: local_coordinate(3)
    integer :: temporary_connectivity(HEX20_NNODE,2)
    integer :: element, local_node, global_node, node_count

    call hex20_natural_coordinates(natural_nodes)
    temporary_coordinates = 0.0_dp
    temporary_connectivity = 0
    node_count = 0
    do element = 1, 2
      do local_node = 1, HEX20_NNODE
        local_coordinate(1) = 0.5_dp*(natural_nodes(1,local_node)+1.0_dp)+real(element-1,dp)
        local_coordinate(2) = 0.5_dp*(natural_nodes(2,local_node)+1.0_dp)
        local_coordinate(3) = 0.5_dp*(natural_nodes(3,local_node)+1.0_dp)
        global_node = locate_coordinate(temporary_coordinates,node_count,local_coordinate)
        if (global_node == 0) then
          node_count = node_count+1
          temporary_coordinates(:,node_count) = local_coordinate
          global_node = node_count
        end if
        temporary_connectivity(local_node,element) = global_node
      end do
    end do
    allocate(mesh%coordinates(3,node_count),mesh%connectivity(HEX20_NNODE,2),mesh%material_id(2))
    mesh%coordinates = temporary_coordinates(:,1:node_count)
    mesh%connectivity = temporary_connectivity
    mesh%material_id = 1
  end subroutine build_two_element_mesh

  integer function locate_coordinate(coordinates,node_count,target) result(node_id)
    real(dp), intent(in) :: coordinates(:,:),target(3)
    integer, intent(in) :: node_count
    integer :: node
    node_id = 0
    do node = 1, node_count
      if (maxval(abs(coordinates(:,node)-target)) < 1.0e-14_dp) then
        node_id = node
        return
      end if
    end do
  end function locate_coordinate

  subroutine define_global_state(mesh,state)
    type(MeshType), intent(in) :: mesh
    type(StateType), intent(inout) :: state
    real(dp) :: x,y,z
    integer :: node
    state%x = 0.0_dp
    do node = 1, mesh%number_of_nodes()
      x=mesh%coordinates(1,node); y=mesh%coordinates(2,node); z=mesh%coordinates(3,node)
      state%x(global_dof(node,F_U1))=1.0e-4_dp*x+0.2e-4_dp*y-0.1e-4_dp*z
      state%x(global_dof(node,F_U2))=-0.3e-4_dp*x+0.8e-4_dp*y+0.2e-4_dp*z
      state%x(global_dof(node,F_U3))=0.1e-4_dp*x-0.2e-4_dp*y+0.6e-4_dp*z
      state%x(global_dof(node,F_PHI))=-1.2e5_dp*x+0.8e5_dp*y-0.6e5_dp*z
      state%x(global_dof(node,F_ETA_N))=0.12_dp+0.05_dp*x-0.03_dp*y+0.02_dp*z
      state%x(global_dof(node,F_ETA_P))=-0.18_dp-0.04_dp*x+0.06_dp*y-0.01_dp*z
      state%x(global_dof(node,F_THETA))=45.0_dp+3.0_dp*x-2.0_dp*y+1.5_dp*z
    end do
  end subroutine define_global_state

  subroutine define_boundary_conditions(mesh,boundary_conditions)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), allocatable, intent(out) :: boundary_conditions(:)
    integer :: left_node,right_node,boundary
    real(dp) :: target(3)
    target=[0.0_dp,0.0_dp,0.0_dp]
    left_node=find_node(mesh,target)
    target=[2.0_dp,0.0_dp,0.0_dp]
    right_node=find_node(mesh,target)
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
    boundary_conditions(8)%value_start=0.0_dp; boundary_conditions(8)%value_end=-2.4e5_dp
  end subroutine define_boundary_conditions

  integer function find_node(mesh,target) result(node_id)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: target(3)
    integer :: node
    node_id=0
    do node=1,mesh%number_of_nodes()
      if(maxval(abs(mesh%coordinates(:,node)-target))<1.0e-14_dp) then
        node_id=node
        return
      end if
    end do
  end function find_node

  subroutine build_global_field_indices(node_count,indices_u,indices_phi,indices_n,indices_p,indices_T)
    integer,intent(in)::node_count
    integer,allocatable,intent(out)::indices_u(:),indices_phi(:),indices_n(:),indices_p(:),indices_T(:)
    integer::node,component,counter
    allocate(indices_u(3*node_count),indices_phi(node_count),indices_n(node_count), &
      indices_p(node_count),indices_T(node_count))
    counter=0
    do node=1,node_count
      do component=F_U1,F_U3
        counter=counter+1; indices_u(counter)=global_dof(node,component)
      end do
      indices_phi(node)=global_dof(node,F_PHI)
      indices_n(node)=global_dof(node,F_ETA_N)
      indices_p(node)=global_dof(node,F_ETA_P)
      indices_T(node)=global_dof(node,F_THETA)
    end do
  end subroutine build_global_field_indices

  real(dp) function check_csr_structure(context) result(error)
    type(SparseAssemblyContextType),intent(in)::context
    integer::row,position
    error=0.0_dp
    if(context%ia(1)/=1 .or. context%ia(context%number_of_equations+1)/= &
      context%number_of_nonzeros+1) error=1.0_dp
    do row=1,context%number_of_equations
      if(find_csr_position(context,row,row)==0) error=1.0_dp
      do position=context%ia(row)+1,context%ia(row+1)-1
        if(context%ja(position)<=context%ja(position-1)) error=1.0_dp
      end do
    end do
  end function check_csr_structure

  real(dp) function check_dirichlet_rows(context,system,fixed) result(error)
    type(SparseAssemblyContextType),intent(in)::context
    type(GlobalSystemType),intent(in)::system
    logical,intent(in)::fixed(:)
    integer::row,position
    real(dp)::expected
    error=0.0_dp
    do row=1,context%number_of_equations
      if(.not.fixed(row)) cycle
      error=max(error,abs(system%residual(row)))
      do position=context%ia(row),context%ia(row+1)-1
        expected=0.0_dp
        if(context%ja(position)==row) expected=1.0_dp
        error=max(error,abs(system%jacobian_values(position)-expected))
      end do
    end do
  end function check_dirichlet_rows

  real(dp) function check_shared_node_scatter(mesh,material,state,load,context,system,status) result(error)
    type(MeshType),intent(in)::mesh
    type(MaterialType),intent(in)::material
    type(StateType),intent(in)::state
    type(LoadStateType),intent(in)::load
    type(SparseAssemblyContextType),intent(in)::context
    type(GlobalSystemType),intent(in)::system
    type(SolverStatusType),intent(out)::status
    type(ElementRequestType)::request
    type(ElementDiagnosticsType)::diagnostics
    type(SolverStatusType)::local_status
    real(dp)::coordinates_e(3,HEX20_NNODE),state_e(HEX20_NDOF),target(3)
    real(dp)::residual_e(HEX20_NDOF),jacobian_e(HEX20_NDOF,HEX20_NDOF),manual
    integer::shared_node,element,local_node
    call clear_status(status)
    target=[1.0_dp,0.0_dp,0.0_dp]
    shared_node=find_node(mesh,target)
    manual=0.0_dp; request%need_jacobian=.false.
    do element=1,mesh%number_of_elements()
      coordinates_e=mesh%coordinates(:,mesh%connectivity(:,element))
      state_e=state%x(context%element_dof(:,element))
      call evaluate_hex20_element(material,coordinates_e,state_e,load,request,residual_e, &
        jacobian_e,diagnostics,local_status)
      if(.not.local_status%is_ok()) then
        call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); error=1.0_dp; return
      end if
      do local_node=1,HEX20_NNODE
        if(mesh%connectivity(local_node,element)==shared_node) &
          manual=manual+residual_e(local_dof(local_node,F_THETA))
      end do
    end do
    error=abs(manual-system%residual(global_dof(shared_node,F_THETA))) / max(1.0_dp,abs(manual))
  end function check_shared_node_scatter

  real(dp) function relative_error(a,b)
    real(dp),intent(in)::a(:),b(:)
    relative_error=vector_norm(a-b)/max(1.0e-30_dp,vector_norm(a),vector_norm(b))
  end function relative_error

  real(dp) function vector_norm(vector)
    real(dp),intent(in)::vector(:)
    vector_norm=sqrt(sum(vector*vector))
  end function vector_norm

  subroutine print_metric(name,value)
    character(len=*),intent(in)::name
    real(dp),intent(in)::value
    write(*,'(a,es24.16e3)')trim(name)//'=',value
  end subroutine print_metric

end module global_assembly_check_m

module reference_bar_problem_m
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE, F_U1, F_U2, F_U3, F_PHI, F_ETA_N, &
    F_ETA_P, F_THETA
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INPUT_ERROR
  use input_schema_m, only : CaseConfigType, PROBLEM_REFERENCE_BAR, &
    PROBLEM_PLAIN_PLATE, PROBLEM_CENTER_HOLE_PLATE, MATERIAL_TPS_G08_REFERENCE
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType
  use equilibrium_m, only : compute_equilibrium_densities
  use reference_material_m, only : define_tps_g08_reference_material
  use dirichlet_m, only : DirichletBCType
  use boundary_topology_m, only : BoundaryTopologyAuditType, audit_boundary_topology
  implicit none
  private

  interface initialize_voltage_driven_problem
    module procedure initialize_reference_bar_problem
  end interface initialize_voltage_driven_problem

  public :: initialize_reference_bar_problem, initialize_voltage_driven_problem

contains

  subroutine initialize_reference_bar_problem(config,mesh,materials,state,start_load, &
                                              boundary_conditions,status)
    type(CaseConfigType), intent(in) :: config
    type(MeshType), intent(in) :: mesh
    type(MaterialType), allocatable, intent(out) :: materials(:)
    type(StateType), intent(out) :: state
    type(LoadStateType), intent(out) :: start_load
    type(DirichletBCType), allocatable, intent(out) :: boundary_conditions(:)
    type(SolverStatusType), intent(out) :: status
    integer, allocatable :: left_nodes(:),right_nodes(:)
    type(SolverStatusType) :: local_status
    type(BoundaryTopologyAuditType) :: topology
    real(dp) :: n0,p0,x_min,x_max,tolerance

    call clear_status(status)
    if ((trim(config%problem_type) /= PROBLEM_REFERENCE_BAR .and. &
         trim(config%problem_type) /= PROBLEM_PLAIN_PLATE .and. &
         trim(config%problem_type) /= PROBLEM_CENTER_HOLE_PLATE) .or. &
         trim(config%material_profile) /= MATERIAL_TPS_G08_REFERENCE) then
      call set_status(status,STATUS_INPUT_ERROR, &
        'Reference-bar setup received an unsupported problem or material profile.')
      return
    end if
    if (mesh%number_of_nodes() < 1 .or. mesh%number_of_elements() < 1) then
      call set_status(status,STATUS_INPUT_ERROR,'Reference-bar mesh is empty.')
      return
    end if
    if (any(mesh%material_id /= 1)) then
      call set_status(status,STATUS_INPUT_ERROR, &
        'TPS_G08_REFERENCE requires material identifier 1 for every element.')
      return
    end if
    call audit_boundary_topology(mesh,config%hole_radius,topology,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,local_status%code,trim(local_status%message)); return
    end if
    if (topology%unknown_face_count /= 0) then
      call set_status(status,STATUS_INPUT_ERROR, &
        'Voltage-driven mesh contains an unclassified exterior boundary face.')
      return
    end if
    if (trim(config%problem_type) == PROBLEM_CENTER_HOLE_PLATE) then
      if (topology%hole_face_count < 1 .or. &
          topology%maximum_hole_radius_error > topology%classification_tolerance) then
        call set_status(status,STATUS_INPUT_ERROR, &
          'CENTER_HOLE_PLATE does not contain the declared cylindrical hole boundary.')
        return
      end if
    else if (topology%hole_face_count /= 0) then
      call set_status(status,STATUS_INPUT_ERROR, &
        'A non-hole problem unexpectedly contains a classified hole boundary.')
      return
    end if

    allocate(materials(1))
    call define_tps_g08_reference_material(materials(1))
    call compute_equilibrium_densities(materials(1)%donor_concentration, &
      materials(1)%acceptor_concentration,materials(1)%intrinsic_concentration, &
      n0,p0,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,local_status%code,trim(local_status%message)); return
    end if
    materials(1)%reference_n=n0
    materials(1)%reference_p=p0

    x_min=minval(mesh%coordinates(1,:))
    x_max=maxval(mesh%coordinates(1,:))
    tolerance=max(1.0e-12_dp*abs(x_max-x_min),1.0e-18_dp)
    if (x_max-x_min <= tolerance) then
      call set_status(status,STATUS_INPUT_ERROR,'Reference-bar x extent is degenerate.')
      return
    end if
    call collect_face_nodes(mesh,x_min,tolerance,left_nodes)
    call collect_face_nodes(mesh,x_max,tolerance,right_nodes)
    if (size(left_nodes) < 8 .or. size(right_nodes) < 8) then
      call set_status(status,STATUS_INPUT_ERROR, &
        'Reference-bar contact faces must contain at least one complete HEX20 face.')
      return
    end if

    allocate(boundary_conditions(11))
    call set_boundary(boundary_conditions(1),F_U1,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(2),F_U2,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(3),F_U3,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(4),F_PHI,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(5),F_PHI,right_nodes,0.0_dp,config%target_voltage)
    call set_boundary(boundary_conditions(6),F_ETA_N,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(7),F_ETA_N,right_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(8),F_ETA_P,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(9),F_ETA_P,right_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(10),F_THETA,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(11),F_THETA,right_nodes,0.0_dp,0.0_dp)

    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()))
    state%x=0.0_dp
    start_load%voltage_factor=0.0_dp
    start_load%joule_factor=0.0_dp
  end subroutine initialize_reference_bar_problem

  subroutine collect_face_nodes(mesh,x_coordinate,tolerance,nodes)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: x_coordinate,tolerance
    integer, allocatable, intent(out) :: nodes(:)
    integer :: node,entry

    allocate(nodes(count(abs(mesh%coordinates(1,:)-x_coordinate) <= tolerance)))
    entry=0
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(1,node)-x_coordinate) > tolerance) cycle
      entry=entry+1
      nodes(entry)=node
    end do
  end subroutine collect_face_nodes

  subroutine set_boundary(boundary,field,nodes,value_start,value_end)
    type(DirichletBCType), intent(out) :: boundary
    integer, intent(in) :: field,nodes(:)
    real(dp), intent(in) :: value_start,value_end

    boundary%field_id=field
    allocate(boundary%node_ids(size(nodes)))
    boundary%node_ids=nodes
    boundary%value_start=value_start
    boundary%value_end=value_end
  end subroutine set_boundary

end module reference_bar_problem_m

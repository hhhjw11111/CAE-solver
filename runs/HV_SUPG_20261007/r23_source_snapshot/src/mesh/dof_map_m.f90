module dof_map_m
  use field_layout_m, only : NDOF_NODE, global_dof, local_dof
  use mesh_types_m, only : HEX20_NNODE
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INVALID_ARGUMENT
  implicit none
  private

  integer, parameter, public :: HEX20_NDOF = HEX20_NNODE * NDOF_NODE
  public :: build_element_dof_map

contains

  subroutine build_element_dof_map(connectivity, element_dof, status)
    integer, intent(in) :: connectivity(HEX20_NNODE)
    integer, intent(out) :: element_dof(HEX20_NDOF)
    type(SolverStatusType), intent(out) :: status
    integer :: node, field

    call clear_status(status)
    if (any(connectivity <= 0)) then
      element_dof = 0
      call set_status(status, STATUS_INVALID_ARGUMENT, &
        'Element connectivity contains a non-positive node identifier.')
      return
    end if

    do node = 1, HEX20_NNODE
      do field = 1, NDOF_NODE
        element_dof(local_dof(node,field)) = global_dof(connectivity(node),field)
      end do
    end do
  end subroutine build_element_dof_map

end module dof_map_m

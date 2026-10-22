module field_layout_m
  implicit none
  private

  integer, parameter, public :: NDOF_NODE = 7
  integer, parameter, public :: F_U1 = 1
  integer, parameter, public :: F_U2 = 2
  integer, parameter, public :: F_U3 = 3
  integer, parameter, public :: F_PHI = 4
  integer, parameter, public :: F_ETA_N = 5
  integer, parameter, public :: F_ETA_P = 6
  integer, parameter, public :: F_THETA = 7

  character(len=12), parameter, public :: FIELD_NAMES(NDOF_NODE) = [ &
    character(len=12) :: 'u1', 'u2', 'u3', 'phi', 'eta_n', 'eta_p', 'theta' ]

  public :: global_dof, local_dof, decode_global_dof

contains

  pure integer function global_dof(node_id, field_id) result(index)
    integer, intent(in) :: node_id, field_id
    index = (node_id - 1) * NDOF_NODE + field_id
  end function global_dof

  pure integer function local_dof(local_node_id, field_id) result(index)
    integer, intent(in) :: local_node_id, field_id
    index = (local_node_id - 1) * NDOF_NODE + field_id
  end function local_dof

  pure subroutine decode_global_dof(index, node_id, field_id)
    integer, intent(in) :: index
    integer, intent(out) :: node_id, field_id
    node_id = (index - 1) / NDOF_NODE + 1
    field_id = mod(index - 1, NDOF_NODE) + 1
  end subroutine decode_global_dof

end module field_layout_m

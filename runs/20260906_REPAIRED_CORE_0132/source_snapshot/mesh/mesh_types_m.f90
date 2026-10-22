module mesh_types_m
  use kinds_m, only : dp
  implicit none
  private

  integer, parameter, public :: HEX20_NNODE = 20
  integer, parameter, public :: HEX20_FACE_NNODE = 8

  type, public :: BoundaryFaceType
    integer :: owner_element = 0
    integer :: local_face = 0
    integer :: boundary_set = 0
    integer :: node_ids(HEX20_FACE_NNODE) = 0
  end type BoundaryFaceType

  type, public :: MeshType
    real(dp), allocatable :: coordinates(:,:)       ! (3, number_of_nodes)
    integer, allocatable :: connectivity(:,:)       ! (20, number_of_elements)
    integer, allocatable :: material_id(:)
    type(BoundaryFaceType), allocatable :: boundary_faces(:)
  contains
    procedure :: number_of_nodes
    procedure :: number_of_elements
  end type MeshType

contains

  integer function number_of_nodes(this)
    class(MeshType), intent(in) :: this
    if (allocated(this%coordinates)) then
      number_of_nodes = size(this%coordinates, 2)
    else
      number_of_nodes = 0
    end if
  end function number_of_nodes

  integer function number_of_elements(this)
    class(MeshType), intent(in) :: this
    if (allocated(this%connectivity)) then
      number_of_elements = size(this%connectivity, 2)
    else
      number_of_elements = 0
    end if
  end function number_of_elements

end module mesh_types_m

module mesh_checks_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INPUT_ERROR, STATUS_GEOMETRY_ERROR
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use hex20_shape_m, only : hex20_geometry
  use hex20_quadrature_m, only : HEX20_NGAUSS, gauss_rule_3x3x3
  implicit none
  private

  public :: validate_mesh

contains

  subroutine validate_mesh(mesh, minimum_determinant, total_volume, status)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(out) :: minimum_determinant, total_volume
    type(SolverStatusType), intent(out) :: status
    real(dp) :: points(3,HEX20_NGAUSS), weights(HEX20_NGAUSS)
    real(dp) :: coordinates_e(3,HEX20_NNODE)
    real(dp) :: shape(HEX20_NNODE), derivative_global(3,HEX20_NNODE), determinant
    type(SolverStatusType) :: local_status
    integer :: element, point, node_a, node_b, node_count, element_count

    call clear_status(status)
    minimum_determinant = huge(1.0_dp)
    total_volume = 0.0_dp

    if (.not. allocated(mesh%coordinates) .or. .not. allocated(mesh%connectivity)) then
      call set_status(status, STATUS_INPUT_ERROR, 'Mesh coordinates or connectivity are not allocated.')
      return
    end if
    if (size(mesh%coordinates,1) /= 3 .or. size(mesh%connectivity,1) /= HEX20_NNODE) then
      call set_status(status, STATUS_INPUT_ERROR, 'Mesh dimensions are not compatible with 3D HEX20.')
      return
    end if

    node_count = mesh%number_of_nodes()
    element_count = mesh%number_of_elements()
    if (node_count < HEX20_NNODE .or. element_count < 1) then
      call set_status(status, STATUS_INPUT_ERROR, 'Mesh must contain at least 20 nodes and one element.')
      return
    end if
    if (.not. all(ieee_is_finite(mesh%coordinates))) then
      call set_status(status, STATUS_INPUT_ERROR, 'Mesh contains a non-finite coordinate.')
      return
    end if
    if (any(mesh%connectivity < 1) .or. any(mesh%connectivity > node_count)) then
      call set_status(status, STATUS_INPUT_ERROR, 'Mesh connectivity references an invalid node.')
      return
    end if

    do element = 1, element_count
      do node_a = 1, HEX20_NNODE-1
        do node_b = node_a+1, HEX20_NNODE
          if (mesh%connectivity(node_a,element) == mesh%connectivity(node_b,element)) then
            call set_status(status, STATUS_INPUT_ERROR, 'HEX20 element contains a repeated node.')
            return
          end if
        end do
      end do
    end do

    call gauss_rule_3x3x3(points, weights)
    do element = 1, element_count
      coordinates_e = mesh%coordinates(:,mesh%connectivity(:,element))
      do point = 1, HEX20_NGAUSS
        call hex20_geometry(coordinates_e, points(1,point), points(2,point), &
          points(3,point), shape, derivative_global, determinant, local_status)
        if (.not. local_status%is_ok()) then
          call set_status(status, STATUS_GEOMETRY_ERROR, &
            'Invalid geometry in element ' // integer_to_string(element) // &
            ', Gauss point ' // integer_to_string(point) // '.')
          return
        end if
        minimum_determinant = min(minimum_determinant, determinant)
        total_volume = total_volume + weights(point) * determinant
      end do
      call validate_element_boundary_trace(coordinates_e,element,local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status,local_status%code,trim(local_status%message))
        return
      end if
    end do
  end subroutine validate_mesh

  subroutine validate_element_boundary_trace(coordinates,element,status)
    real(dp), intent(in) :: coordinates(3,HEX20_NNODE)
    integer, intent(in) :: element
    type(SolverStatusType), intent(out) :: status
    real(dp) :: points(3,HEX20_NGAUSS),weights(HEX20_NGAUSS)
    real(dp) :: shape(HEX20_NNODE),derivative_global(3,HEX20_NNODE)
    real(dp) :: determinant
    integer :: face,axis,point
    real(dp) :: rst(3)

    call clear_status(status)
    call gauss_rule_3x3x3(points,weights)
    do face=1,6
      axis=(face+1)/2
      do point=1,HEX20_NGAUSS
        rst=points(:,point)
        rst(axis)=real(2*mod(face+1,2)-1,dp)
        call hex20_geometry(coordinates,rst(1),rst(2),rst(3),shape, &
          derivative_global,determinant,status)
        if (.not. status%is_ok()) then
          call set_status(status,STATUS_GEOMETRY_ERROR, &
            'Invalid HEX20 boundary-trace geometry in element '// &
            trim(integer_to_string(element))//', face '// &
            trim(integer_to_string(face))//', boundary point '// &
            trim(integer_to_string(point))//'.')
          return
        end if
      end do
    end do
  end subroutine validate_element_boundary_trace

  function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

end module mesh_checks_m

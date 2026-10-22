module hex20_shape_m
  use kinds_m, only : dp
  use constants_m, only : GEOMETRY_TOL
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_GEOMETRY_ERROR
  use mesh_types_m, only : HEX20_NNODE
  implicit none
  private

  real(dp), parameter :: NODE_R(HEX20_NNODE) = [ &
    -1.0_dp,  1.0_dp,  1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp,  1.0_dp, -1.0_dp, &
     0.0_dp,  1.0_dp,  0.0_dp, -1.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, -1.0_dp, &
    -1.0_dp,  1.0_dp,  1.0_dp, -1.0_dp ]
  real(dp), parameter :: NODE_S(HEX20_NNODE) = [ &
    -1.0_dp, -1.0_dp,  1.0_dp,  1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp,  1.0_dp, &
    -1.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, -1.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, &
    -1.0_dp, -1.0_dp,  1.0_dp,  1.0_dp ]
  real(dp), parameter :: NODE_T(HEX20_NNODE) = [ &
    -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp,  1.0_dp,  1.0_dp,  1.0_dp, &
    -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp,  1.0_dp,  1.0_dp,  1.0_dp,  1.0_dp, &
     0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp ]

  public :: hex20_shape, hex20_geometry, hex20_natural_coordinates

contains

  subroutine hex20_natural_coordinates(coordinates)
    real(dp), intent(out) :: coordinates(3,HEX20_NNODE)
    coordinates(1,:) = NODE_R
    coordinates(2,:) = NODE_S
    coordinates(3,:) = NODE_T
  end subroutine hex20_natural_coordinates

  subroutine hex20_shape(r, s, t, shape, derivative_natural)
    real(dp), intent(in) :: r, s, t
    real(dp), intent(out) :: shape(HEX20_NNODE)
    real(dp), intent(out) :: derivative_natural(3,HEX20_NNODE)
    real(dp) :: a, b, c, ar, bs, ct, linear_factor
    integer :: node

    shape = 0.0_dp
    derivative_natural = 0.0_dp

    do node = 1, 8
      a = NODE_R(node)
      b = NODE_S(node)
      c = NODE_T(node)
      ar = 1.0_dp + a*r
      bs = 1.0_dp + b*s
      ct = 1.0_dp + c*t
      linear_factor = a*r + b*s + c*t - 2.0_dp
      shape(node) = 0.125_dp * ar * bs * ct * linear_factor
      derivative_natural(1,node) = 0.125_dp * a * bs * ct * (linear_factor + ar)
      derivative_natural(2,node) = 0.125_dp * b * ar * ct * (linear_factor + bs)
      derivative_natural(3,node) = 0.125_dp * c * ar * bs * (linear_factor + ct)
    end do

    do node = 9, HEX20_NNODE
      a = NODE_R(node)
      b = NODE_S(node)
      c = NODE_T(node)
      if (a == 0.0_dp) then
        shape(node) = 0.25_dp * (1.0_dp-r*r) * (1.0_dp+b*s) * (1.0_dp+c*t)
        derivative_natural(1,node) = -0.5_dp*r * (1.0_dp+b*s) * (1.0_dp+c*t)
        derivative_natural(2,node) = 0.25_dp*(1.0_dp-r*r)*b*(1.0_dp+c*t)
        derivative_natural(3,node) = 0.25_dp*(1.0_dp-r*r)*c*(1.0_dp+b*s)
      else if (b == 0.0_dp) then
        shape(node) = 0.25_dp * (1.0_dp-s*s) * (1.0_dp+a*r) * (1.0_dp+c*t)
        derivative_natural(1,node) = 0.25_dp*(1.0_dp-s*s)*a*(1.0_dp+c*t)
        derivative_natural(2,node) = -0.5_dp*s * (1.0_dp+a*r) * (1.0_dp+c*t)
        derivative_natural(3,node) = 0.25_dp*(1.0_dp-s*s)*c*(1.0_dp+a*r)
      else
        shape(node) = 0.25_dp * (1.0_dp-t*t) * (1.0_dp+a*r) * (1.0_dp+b*s)
        derivative_natural(1,node) = 0.25_dp*(1.0_dp-t*t)*a*(1.0_dp+b*s)
        derivative_natural(2,node) = 0.25_dp*(1.0_dp-t*t)*b*(1.0_dp+a*r)
        derivative_natural(3,node) = -0.5_dp*t * (1.0_dp+a*r) * (1.0_dp+b*s)
      end if
    end do
  end subroutine hex20_shape

  subroutine hex20_geometry(coordinates, r, s, t, shape, derivative_global, &
                            determinant, status)
    real(dp), intent(in) :: coordinates(3,HEX20_NNODE)
    real(dp), intent(in) :: r, s, t
    real(dp), intent(out) :: shape(HEX20_NNODE)
    real(dp), intent(out) :: derivative_global(3,HEX20_NNODE)
    real(dp), intent(out) :: determinant
    type(SolverStatusType), intent(out) :: status
    real(dp) :: derivative_natural(3,HEX20_NNODE)
    real(dp) :: jacobian(3,3), inverse_jacobian(3,3)
    real(dp) :: jacobian_scale, determinant_tolerance

    call clear_status(status)
    call hex20_shape(r, s, t, shape, derivative_natural)
    jacobian = matmul(derivative_natural, transpose(coordinates))
    call inverse_3x3(jacobian, inverse_jacobian, determinant)
    jacobian_scale = maxval(abs(jacobian))
    determinant_tolerance = GEOMETRY_TOL * jacobian_scale**3
    if (jacobian_scale <= tiny(1.0_dp) .or. determinant <= determinant_tolerance) then
      derivative_global = 0.0_dp
      call set_status(status, STATUS_GEOMETRY_ERROR, &
        'HEX20 element has a non-positive or relatively near-zero Jacobian determinant.')
      return
    end if
    derivative_global = matmul(inverse_jacobian, derivative_natural)
  end subroutine hex20_geometry

  pure subroutine inverse_3x3(matrix, inverse_matrix, determinant)
    real(dp), intent(in) :: matrix(3,3)
    real(dp), intent(out) :: inverse_matrix(3,3)
    real(dp), intent(out) :: determinant

    determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
                - matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
                + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
    if (abs(determinant) <= tiny(1.0_dp)) then
      inverse_matrix = 0.0_dp
      return
    end if

    inverse_matrix(1,1) =  (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant
    inverse_matrix(1,2) = -(matrix(1,2)*matrix(3,3)-matrix(1,3)*matrix(3,2))/determinant
    inverse_matrix(1,3) =  (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
    inverse_matrix(2,1) = -(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))/determinant
    inverse_matrix(2,2) =  (matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant
    inverse_matrix(2,3) = -(matrix(1,1)*matrix(2,3)-matrix(1,3)*matrix(2,1))/determinant
    inverse_matrix(3,1) =  (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant
    inverse_matrix(3,2) = -(matrix(1,1)*matrix(3,2)-matrix(1,2)*matrix(3,1))/determinant
    inverse_matrix(3,3) =  (matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
  end subroutine inverse_3x3

end module hex20_shape_m

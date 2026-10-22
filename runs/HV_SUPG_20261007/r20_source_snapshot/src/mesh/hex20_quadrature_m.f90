module hex20_quadrature_m
  use kinds_m, only : dp
  implicit none
  private

  integer, parameter, public :: HEX20_NGAUSS = 27
  public :: gauss_rule_3x3x3

contains

  subroutine gauss_rule_3x3x3(points, weights)
    real(dp), intent(out) :: points(3,HEX20_NGAUSS)
    real(dp), intent(out) :: weights(HEX20_NGAUSS)
    real(dp) :: location(3), weight_1d(3)
    integer :: ir, is, it, point_id

    location = [ -sqrt(3.0_dp/5.0_dp), 0.0_dp, sqrt(3.0_dp/5.0_dp) ]
    weight_1d = [ 5.0_dp/9.0_dp, 8.0_dp/9.0_dp, 5.0_dp/9.0_dp ]

    point_id = 0
    do ir = 1, 3
      do is = 1, 3
        do it = 1, 3
          point_id = point_id + 1
          points(1,point_id) = location(ir)
          points(2,point_id) = location(is)
          points(3,point_id) = location(it)
          weights(point_id) = weight_1d(ir) * weight_1d(is) * weight_1d(it)
        end do
      end do
    end do
  end subroutine gauss_rule_3x3x3

end module hex20_quadrature_m

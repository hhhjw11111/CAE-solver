module hex20_quadrature_m
  use kinds_m, only : dp
  implicit none
  private

  integer, parameter, public :: HEX20_NGAUSS = 27
  public :: gauss_rule_3x3x3, gauss_rule_5x5x5

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

  subroutine gauss_rule_5x5x5(points,weights)
    real(dp), intent(out) :: points(3,125),weights(125)
    real(dp) :: x(5),w(5),inner,outer
    integer :: i,j,k,n
    inner=sqrt(5.0_dp-2.0_dp*sqrt(10.0_dp/7.0_dp))/3.0_dp
    outer=sqrt(5.0_dp+2.0_dp*sqrt(10.0_dp/7.0_dp))/3.0_dp
    x=[-outer,-inner,0.0_dp,inner,outer]
    w=[(322.0_dp-13.0_dp*sqrt(70.0_dp))/900.0_dp, &
       (322.0_dp+13.0_dp*sqrt(70.0_dp))/900.0_dp,128.0_dp/225.0_dp, &
       (322.0_dp+13.0_dp*sqrt(70.0_dp))/900.0_dp, &
       (322.0_dp-13.0_dp*sqrt(70.0_dp))/900.0_dp]
    n=0
    do i=1,5
      do j=1,5
        do k=1,5
          n=n+1
          points(1,n)=x(i)
          points(2,n)=x(j)
          points(3,n)=x(k)
          weights(n)=w(i)*w(j)*w(k)
        end do
      end do
    end do
  end subroutine gauss_rule_5x5x5
end module hex20_quadrature_m

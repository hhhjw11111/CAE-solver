module recombination_m
  use kinds_m, only : dp
  implicit none
  private

  public :: evaluate_recombination_none

contains

  pure subroutine evaluate_recombination_none(n, p, temperature, rate, rate_n, rate_p, rate_T)
    real(dp), intent(in) :: n, p, temperature
    real(dp), intent(out) :: rate, rate_n, rate_p, rate_T

    rate = 0.0_dp*n + 0.0_dp*p + 0.0_dp*temperature
    rate_n = 0.0_dp
    rate_p = 0.0_dp
    rate_T = 0.0_dp
  end subroutine evaluate_recombination_none

end module recombination_m

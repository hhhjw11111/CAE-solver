module equilibrium_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INVALID_ARGUMENT
  implicit none
  private

  public :: compute_equilibrium_densities

contains

  subroutine compute_equilibrium_densities(donor, acceptor, intrinsic, n0, p0, status)
    real(dp), intent(in) :: donor, acceptor, intrinsic
    real(dp), intent(out) :: n0, p0
    type(SolverStatusType), intent(out) :: status
    real(dp) :: net_doping, discriminant

    call clear_status(status)
    n0 = 0.0_dp
    p0 = 0.0_dp
    if (donor < 0.0_dp .or. acceptor < 0.0_dp .or. intrinsic <= 0.0_dp) then
      call set_status(status, STATUS_INVALID_ARGUMENT, &
        'Doping must be non-negative and intrinsic concentration must be positive.')
      return
    end if

    net_doping = donor - acceptor
    discriminant = sqrt(net_doping*net_doping + 4.0_dp*intrinsic*intrinsic)
    if (net_doping >= 0.0_dp) then
      n0 = 0.5_dp * (net_doping + discriminant)
      p0 = intrinsic*intrinsic / n0
    else
      p0 = 0.5_dp * (-net_doping + discriminant)
      n0 = intrinsic*intrinsic / p0
    end if
  end subroutine compute_equilibrium_densities

end module equilibrium_m

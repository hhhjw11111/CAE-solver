module transport_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE, BOLTZMANN_CONSTANT
  use model_types_m, only : MaterialType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INVALID_ARGUMENT
  implicit none
  private

  public :: evaluate_transport

contains

  subroutine evaluate_transport(material, temperature, mobility_n, mobility_p, &
                                diffusion_n, diffusion_p, mobility_n_T, mobility_p_T, &
                                diffusion_n_T, diffusion_p_T, status)
    type(MaterialType), intent(in) :: material
    real(dp), intent(in) :: temperature
    real(dp), intent(out) :: mobility_n(3,3), mobility_p(3,3)
    real(dp), intent(out) :: diffusion_n(3,3), diffusion_p(3,3)
    real(dp), intent(out) :: mobility_n_T(3,3), mobility_p_T(3,3)
    real(dp), intent(out) :: diffusion_n_T(3,3), diffusion_p_T(3,3)
    type(SolverStatusType), intent(out) :: status
    real(dp) :: temperature_ratio, factor_n, factor_p

    call clear_status(status)
    mobility_n = 0.0_dp
    mobility_p = 0.0_dp
    diffusion_n = 0.0_dp
    diffusion_p = 0.0_dp
    mobility_n_T = 0.0_dp
    mobility_p_T = 0.0_dp
    diffusion_n_T = 0.0_dp
    diffusion_p_T = 0.0_dp

    if (temperature <= 0.0_dp .or. material%reference_temperature <= 0.0_dp) then
      call set_status(status, STATUS_INVALID_ARGUMENT, &
        'Temperature and reference temperature must be positive.')
      return
    end if

    if (material%frozen_transport) then
      mobility_n=material%mobility_n0;mobility_p=material%mobility_p0
      diffusion_n=(BOLTZMANN_CONSTANT*material%reference_temperature/ELEMENTARY_CHARGE)*mobility_n
      diffusion_p=(BOLTZMANN_CONSTANT*material%reference_temperature/ELEMENTARY_CHARGE)*mobility_p
      return ! All four T derivatives were initialized to zero above.
    end if
    temperature_ratio = temperature / material%reference_temperature
    factor_n = temperature_ratio**(-material%beta_n)
    factor_p = temperature_ratio**(-material%beta_p)
    mobility_n = material%mobility_n0 * factor_n
    mobility_p = material%mobility_p0 * factor_p
    mobility_n_T = -(material%beta_n/temperature) * mobility_n
    mobility_p_T = -(material%beta_p/temperature) * mobility_p

    diffusion_n = (BOLTZMANN_CONSTANT*temperature/ELEMENTARY_CHARGE) * mobility_n
    diffusion_p = (BOLTZMANN_CONSTANT*temperature/ELEMENTARY_CHARGE) * mobility_p
    diffusion_n_T = ((1.0_dp-material%beta_n)/temperature) * diffusion_n
    diffusion_p_T = ((1.0_dp-material%beta_p)/temperature) * diffusion_p
  end subroutine evaluate_transport

end module transport_m

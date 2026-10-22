module material_point_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use model_types_m, only : MaterialType, LoadStateType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INVALID_ARGUMENT
  use transport_m, only : evaluate_transport
  use recombination_m, only : evaluate_recombination_none
  use thermo_piezo_m, only : evaluate_thermo_piezo
  implicit none
  private

  integer, parameter, public :: MATERIAL_INPUT_SIZE = 21
  integer, parameter, public :: MATERIAL_RESPONSE_SIZE = 21

  type, public :: MaterialPointStateType
    real(dp) :: strain(6) = 0.0_dp
    real(dp) :: electric_field(3) = 0.0_dp
    real(dp) :: temperature = 300.0_dp
    real(dp) :: temperature_gradient(3) = 0.0_dp
    real(dp) :: eta_n = 0.0_dp
    real(dp) :: eta_n_gradient(3) = 0.0_dp
    real(dp) :: eta_p = 0.0_dp
    real(dp) :: eta_p_gradient(3) = 0.0_dp
  end type MaterialPointStateType

  type, public :: MaterialPointResponseType
    real(dp) :: stress(6) = 0.0_dp
    real(dp) :: electric_displacement(3) = 0.0_dp
    real(dp) :: charge_density = 0.0_dp
    real(dp) :: current_n(3) = 0.0_dp
    real(dp) :: current_p(3) = 0.0_dp
    real(dp) :: heat_flux(3) = 0.0_dp
    real(dp) :: joule_heat = 0.0_dp
    real(dp) :: joule_heat_n = 0.0_dp
    real(dp) :: joule_heat_p = 0.0_dp
    real(dp) :: recombination = 0.0_dp
    real(dp) :: n = 0.0_dp
    real(dp) :: p = 0.0_dp
    real(dp) :: mobility_n(3,3) = 0.0_dp
    real(dp) :: mobility_p(3,3) = 0.0_dp
    real(dp) :: diffusion_n(3,3) = 0.0_dp
    real(dp) :: diffusion_p(3,3) = 0.0_dp
  end type MaterialPointResponseType

  type, public :: MaterialPointTangentType
    real(dp) :: stress_strain(6,6) = 0.0_dp
    real(dp) :: stress_E(6,3) = 0.0_dp
    real(dp) :: stress_T(6) = 0.0_dp
    real(dp) :: displacement_strain(3,6) = 0.0_dp
    real(dp) :: displacement_E(3,3) = 0.0_dp
    real(dp) :: displacement_T(3) = 0.0_dp
    real(dp) :: charge_eta_n = 0.0_dp
    real(dp) :: charge_eta_p = 0.0_dp
    real(dp) :: current_n_E(3,3) = 0.0_dp
    real(dp) :: current_n_eta_n(3) = 0.0_dp
    real(dp) :: current_n_grad_eta_n(3,3) = 0.0_dp
    real(dp) :: current_n_T(3) = 0.0_dp
    real(dp) :: current_p_E(3,3) = 0.0_dp
    real(dp) :: current_p_eta_p(3) = 0.0_dp
    real(dp) :: current_p_grad_eta_p(3,3) = 0.0_dp
    real(dp) :: current_p_T(3) = 0.0_dp
    real(dp) :: heat_flux_grad_T(3,3) = 0.0_dp
    real(dp) :: joule_E(3) = 0.0_dp
    real(dp) :: joule_eta_n = 0.0_dp
    real(dp) :: joule_eta_p = 0.0_dp
    real(dp) :: joule_grad_eta_n(3) = 0.0_dp
    real(dp) :: joule_grad_eta_p(3) = 0.0_dp
    real(dp) :: joule_T = 0.0_dp
    real(dp) :: recombination_eta_n = 0.0_dp
    real(dp) :: recombination_eta_p = 0.0_dp
    real(dp) :: recombination_T = 0.0_dp
  end type MaterialPointTangentType

  public :: evaluate_material_point, pack_material_response, material_tangent_matrix

contains

  subroutine evaluate_material_point(material, state, load, response, tangent, status)
    type(MaterialType), intent(in) :: material
    type(MaterialPointStateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(MaterialPointResponseType), intent(out) :: response
    type(MaterialPointTangentType), intent(out) :: tangent
    type(SolverStatusType), intent(out) :: status
    real(dp) :: mobility_n_T(3,3), mobility_p_T(3,3)
    real(dp) :: diffusion_n_T(3,3), diffusion_p_T(3,3)
    real(dp) :: drift_diffusion_n(3), drift_diffusion_p(3)
    real(dp) :: mobility_field_n(3), mobility_field_p(3)
    real(dp) :: g_n, g_p, rate_n, rate_p, rate_T
    real(dp) :: work_n, work_p

    call clear_status(status)
    response = MaterialPointResponseType()
    tangent = MaterialPointTangentType()
    if (material%reference_n <= 0.0_dp .or. material%reference_p <= 0.0_dp) then
      call set_status(status, STATUS_INVALID_ARGUMENT, &
        'Log-density reference concentrations must be positive.')
      return
    end if
    if (state%temperature <= 0.0_dp) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Material-point temperature must be positive.')
      return
    end if
    if (load%joule_factor < 0.0_dp .or. load%joule_factor > 1.0_dp) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'joule_factor must lie in [0,1].')
      return
    end if

    response%n = material%reference_n * exp(state%eta_n)
    response%p = material%reference_p * exp(state%eta_p)
    if (material%linear_carriers) then
      if (.not. material%frozen_transport) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'LQ requires frozen mobility AND diffusion.');return
      end if
      response%n=material%reference_n*(1.0_dp+state%eta_n)
      response%p=material%reference_p*(1.0_dp+state%eta_p)
      if (response%n<=0.0_dp .or. response%p<=0.0_dp) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'LQ relative density is non-positive.');return
      end if
    end if
    if (.not. ieee_is_finite(response%n) .or. .not. ieee_is_finite(response%p)) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Log-density exponentiation produced a non-finite value.')
      return
    end if

    call evaluate_transport(material, state%temperature, response%mobility_n, &
      response%mobility_p, response%diffusion_n, response%diffusion_p, &
      mobility_n_T, mobility_p_T, diffusion_n_T, diffusion_p_T, status)
    if (.not. status%is_ok()) return

    call evaluate_thermo_piezo(material, state%strain, state%electric_field, &
      state%temperature, response%stress, response%electric_displacement)
    response%charge_density = ELEMENTARY_CHARGE * (response%p-response%n &
      + material%donor_concentration-material%acceptor_concentration)

    mobility_field_n = matmul(response%mobility_n, state%electric_field)
    mobility_field_p = matmul(response%mobility_p, state%electric_field)
    drift_diffusion_n = mobility_field_n &
      + matmul(response%diffusion_n, state%eta_n_gradient)
    drift_diffusion_p = mobility_field_p &
      - matmul(response%diffusion_p, state%eta_p_gradient)
    response%current_n = ELEMENTARY_CHARGE * response%n * drift_diffusion_n
    response%current_p = ELEMENTARY_CHARGE * response%p * drift_diffusion_p
    response%heat_flux = -matmul(material%thermal_conductivity, state%temperature_gradient)

    g_n = dot_product(state%electric_field, mobility_field_n)
    g_p = dot_product(state%electric_field, mobility_field_p)
    response%joule_heat_n = load%joule_factor * ELEMENTARY_CHARGE * response%n * g_n
    response%joule_heat_p = load%joule_factor * ELEMENTARY_CHARGE * response%p * g_p
    response%joule_heat = response%joule_heat_n + response%joule_heat_p
    call evaluate_recombination_none(response%n, response%p, state%temperature, &
      response%recombination, rate_n, rate_p, rate_T)

    tangent%stress_strain = material%elastic
    tangent%stress_E = -transpose(material%piezoelectric)
    tangent%stress_T = -material%thermal_stress
    tangent%displacement_strain = material%piezoelectric
    tangent%displacement_E = material%permittivity
    tangent%displacement_T = material%pyroelectric
    tangent%charge_eta_n = -ELEMENTARY_CHARGE * response%n
    tangent%charge_eta_p =  ELEMENTARY_CHARGE * response%p

    tangent%current_n_E = ELEMENTARY_CHARGE * response%n * response%mobility_n
    tangent%current_n_eta_n = response%current_n
    tangent%current_n_grad_eta_n = ELEMENTARY_CHARGE * response%n * response%diffusion_n
    tangent%current_n_T = ELEMENTARY_CHARGE * response%n * ( &
      matmul(mobility_n_T,state%electric_field) &
      + matmul(diffusion_n_T,state%eta_n_gradient))
    tangent%current_p_E = ELEMENTARY_CHARGE * response%p * response%mobility_p
    tangent%current_p_eta_p = response%current_p
    tangent%current_p_grad_eta_p = -ELEMENTARY_CHARGE * response%p * response%diffusion_p
    tangent%current_p_T = ELEMENTARY_CHARGE * response%p * ( &
      matmul(mobility_p_T,state%electric_field) &
      - matmul(diffusion_p_T,state%eta_p_gradient))

    tangent%heat_flux_grad_T = -material%thermal_conductivity
    tangent%joule_E = load%joule_factor * 2.0_dp * ELEMENTARY_CHARGE * ( &
      response%n*mobility_field_n + response%p*mobility_field_p)
    tangent%joule_eta_n = response%joule_heat_n
    tangent%joule_eta_p = response%joule_heat_p
    tangent%joule_T = load%joule_factor * ELEMENTARY_CHARGE * ( &
      response%n*dot_product(state%electric_field,matmul(mobility_n_T,state%electric_field)) &
      + response%p*dot_product(state%electric_field,matmul(mobility_p_T,state%electric_field)))
    if (material%heat_source_current_work) then
      ! Diagnostic electrical work of the complete drift-diffusion current.
      ! This branch is deliberately not advertised as a thermoelectric closure.
      work_n=dot_product(state%electric_field,response%current_n)
      work_p=dot_product(state%electric_field,response%current_p)
      response%joule_heat_n=load%joule_factor*work_n
      response%joule_heat_p=load%joule_factor*work_p
      response%joule_heat=response%joule_heat_n+response%joule_heat_p
      tangent%joule_E=load%joule_factor*(response%current_n+response%current_p + &
        matmul(transpose(tangent%current_n_E)+transpose(tangent%current_p_E),state%electric_field))
      tangent%joule_eta_n=load%joule_factor*dot_product(state%electric_field,tangent%current_n_eta_n)
      tangent%joule_eta_p=load%joule_factor*dot_product(state%electric_field,tangent%current_p_eta_p)
      tangent%joule_grad_eta_n=load%joule_factor*matmul(transpose(tangent%current_n_grad_eta_n), &
        state%electric_field)
      tangent%joule_grad_eta_p=load%joule_factor*matmul(transpose(tangent%current_p_grad_eta_p), &
        state%electric_field)
      tangent%joule_T=load%joule_factor*dot_product(state%electric_field, &
        tangent%current_n_T+tangent%current_p_T)
    end if
    tangent%recombination_eta_n = rate_n * response%n
    tangent%recombination_eta_p = rate_p * response%p
    tangent%recombination_T = rate_T

    if (material%linear_carriers) then
      ! Exact first-order transport; no eta*E or eta*grad(eta) products.
      response%current_n=ELEMENTARY_CHARGE*material%reference_n*drift_diffusion_n
      response%current_p=ELEMENTARY_CHARGE*material%reference_p*drift_diffusion_p
      tangent%charge_eta_n=-ELEMENTARY_CHARGE*material%reference_n
      tangent%charge_eta_p= ELEMENTARY_CHARGE*material%reference_p
      tangent%current_n_E=ELEMENTARY_CHARGE*material%reference_n*response%mobility_n
      tangent%current_p_E=ELEMENTARY_CHARGE*material%reference_p*response%mobility_p
      tangent%current_n_grad_eta_n=ELEMENTARY_CHARGE*material%reference_n*response%diffusion_n
      tangent%current_p_grad_eta_p=-ELEMENTARY_CHARGE*material%reference_p*response%diffusion_p
      tangent%current_n_eta_n=0.0_dp;tangent%current_p_eta_p=0.0_dp
      tangent%current_n_T=0.0_dp;tangent%current_p_T=0.0_dp
      ! Retained second-order heat is based on reference densities only.
      response%joule_heat_n=load%joule_factor*ELEMENTARY_CHARGE*material%reference_n*g_n
      response%joule_heat_p=load%joule_factor*ELEMENTARY_CHARGE*material%reference_p*g_p
      response%joule_heat=response%joule_heat_n+response%joule_heat_p
      tangent%joule_E=load%joule_factor*2.0_dp*ELEMENTARY_CHARGE*( &
        material%reference_n*mobility_field_n+material%reference_p*mobility_field_p)
      tangent%joule_eta_n=0.0_dp;tangent%joule_eta_p=0.0_dp;tangent%joule_T=0.0_dp
      tangent%joule_grad_eta_n=0.0_dp;tangent%joule_grad_eta_p=0.0_dp
      if (material%heat_source_current_work) then
        response%joule_heat_n=load%joule_factor*dot_product(state%electric_field,response%current_n)
        response%joule_heat_p=load%joule_factor*dot_product(state%electric_field,response%current_p)
        response%joule_heat=response%joule_heat_n+response%joule_heat_p
        tangent%joule_E=load%joule_factor*(response%current_n+response%current_p + &
          matmul(transpose(tangent%current_n_E)+transpose(tangent%current_p_E),state%electric_field))
        tangent%joule_grad_eta_n=load%joule_factor*matmul(transpose(tangent%current_n_grad_eta_n), &
          state%electric_field)
        tangent%joule_grad_eta_p=load%joule_factor*matmul(transpose(tangent%current_p_grad_eta_p), &
          state%electric_field)
      end if
    end if

    if (.not. response_is_finite(response)) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Material-point response contains NaN or Inf.')
    end if
  end subroutine evaluate_material_point

  logical function response_is_finite(response)
    type(MaterialPointResponseType), intent(in) :: response
    response_is_finite = all(ieee_is_finite(response%stress)) .and. &
      all(ieee_is_finite(response%electric_displacement)) .and. &
      ieee_is_finite(response%charge_density) .and. &
      all(ieee_is_finite(response%current_n)) .and. &
      all(ieee_is_finite(response%current_p)) .and. &
      all(ieee_is_finite(response%heat_flux)) .and. &
      ieee_is_finite(response%joule_heat)
  end function response_is_finite

  pure subroutine pack_material_response(response, packed)
    type(MaterialPointResponseType), intent(in) :: response
    real(dp), intent(out) :: packed(MATERIAL_RESPONSE_SIZE)

    packed(1:6) = response%stress
    packed(7:9) = response%electric_displacement
    packed(10) = response%charge_density
    packed(11:13) = response%current_n
    packed(14:16) = response%current_p
    packed(17:19) = response%heat_flux
    packed(20) = response%joule_heat
    packed(21) = response%recombination
  end subroutine pack_material_response

  pure subroutine material_tangent_matrix(tangent, matrix)
    type(MaterialPointTangentType), intent(in) :: tangent
    real(dp), intent(out) :: matrix(MATERIAL_RESPONSE_SIZE,MATERIAL_INPUT_SIZE)

    matrix = 0.0_dp
    matrix(1:6,1:6) = tangent%stress_strain
    matrix(1:6,7:9) = tangent%stress_E
    matrix(1:6,10) = tangent%stress_T
    matrix(7:9,1:6) = tangent%displacement_strain
    matrix(7:9,7:9) = tangent%displacement_E
    matrix(7:9,10) = tangent%displacement_T
    matrix(10,14) = tangent%charge_eta_n
    matrix(10,18) = tangent%charge_eta_p
    matrix(11:13,7:9) = tangent%current_n_E
    matrix(11:13,10) = tangent%current_n_T
    matrix(11:13,14) = tangent%current_n_eta_n
    matrix(11:13,15:17) = tangent%current_n_grad_eta_n
    matrix(14:16,7:9) = tangent%current_p_E
    matrix(14:16,10) = tangent%current_p_T
    matrix(14:16,18) = tangent%current_p_eta_p
    matrix(14:16,19:21) = tangent%current_p_grad_eta_p
    matrix(17:19,11:13) = tangent%heat_flux_grad_T
    matrix(20,7:9) = tangent%joule_E
    matrix(20,10) = tangent%joule_T
    matrix(20,14) = tangent%joule_eta_n
    matrix(20,15:17) = tangent%joule_grad_eta_n
    matrix(20,18) = tangent%joule_eta_p
    matrix(20,19:21) = tangent%joule_grad_eta_p
    matrix(21,10) = tangent%recombination_T
    matrix(21,14) = tangent%recombination_eta_n
    matrix(21,18) = tangent%recombination_eta_p
  end subroutine material_tangent_matrix

end module material_point_m

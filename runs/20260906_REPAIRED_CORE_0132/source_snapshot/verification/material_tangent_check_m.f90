module material_tangent_check_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE, BOLTZMANN_CONSTANT
  use model_types_m, only : MaterialType, LoadStateType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use equilibrium_m, only : compute_equilibrium_densities
  use reference_material_m, only : define_tps_g08_reference_material
  use material_point_m, only : MaterialPointStateType, MaterialPointResponseType, &
    MaterialPointTangentType, evaluate_material_point, pack_material_response, &
    material_tangent_matrix
  implicit none
  private

  integer, parameter :: INPUT_SIZE = 21
  integer, parameter :: RESPONSE_SIZE = 21
  integer, parameter :: NUMBER_OF_GROUPS = 8

  public :: verify_material_point_g02, define_test_material

contains

  subroutine verify_material_point_g02(status)
    type(SolverStatusType), intent(out) :: status
    type(MaterialType) :: material
    type(LoadStateType) :: load
    type(MaterialPointStateType) :: base_state, plus_state, minus_state
    type(MaterialPointResponseType) :: base_response, plus_response, minus_response
    type(MaterialPointTangentType) :: tangent
    type(SolverStatusType) :: local_status
    real(dp) :: analytic(RESPONSE_SIZE,INPUT_SIZE), direction(INPUT_SIZE), raw(INPUT_SIZE)
    real(dp) :: response_plus(RESPONSE_SIZE), response_minus(RESPONSE_SIZE)
    real(dp) :: finite_difference(RESPONSE_SIZE), analytic_action(RESPONSE_SIZE)
    real(dp) :: minimum_error(NUMBER_OF_GROUPS), current_error(NUMBER_OF_GROUPS)
    real(dp) :: step_values(9), step, n0, p0, neutrality_error, mass_action_error
    real(dp) :: einstein_error_n, einstein_error_p
    integer :: index, step_index
    logical :: passed

    call clear_status(status)
    call define_test_material(material)
    call define_test_state(base_state)
    load%voltage_factor = 1.0_dp
    load%joule_factor = 0.73_dp

    call compute_equilibrium_densities(material%donor_concentration, &
      material%acceptor_concentration, material%intrinsic_concentration, &
      n0, p0, local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
      return
    end if
    neutrality_error = abs((p0-n0+material%donor_concentration &
      - material%acceptor_concentration) / max(n0,p0))
    mass_action_error = abs((n0*p0-material%intrinsic_concentration**2) &
      / material%intrinsic_concentration**2)

    call evaluate_material_point(material, base_state, load, base_response, tangent, local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
      return
    end if
    call material_tangent_matrix(tangent, analytic)

    do index = 1, INPUT_SIZE
      raw(index) = sin(0.37_dp*real(index,dp)) + cos(0.19_dp*real(index,dp))
    end do
    raw = raw / sqrt(sum(raw*raw))
    direction = raw * input_scales()
    analytic_action = matmul(analytic, direction)

    step_values = [1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp, 3.0e-4_dp, 1.0e-4_dp, &
                   3.0e-5_dp, 1.0e-5_dp, 3.0e-6_dp, 1.0e-6_dp]
    minimum_error = huge(1.0_dp)
    do step_index = 1, size(step_values)
      step = step_values(step_index)
      call perturb_state(base_state, direction, step, plus_state)
      call perturb_state(base_state, direction, -step, minus_state)
      call evaluate_material_point(material, plus_state, load, plus_response, tangent, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
        return
      end if
      call evaluate_material_point(material, minus_state, load, minus_response, tangent, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
        return
      end if
      call pack_material_response(plus_response, response_plus)
      call pack_material_response(minus_response, response_minus)
      finite_difference = (response_plus-response_minus)/(2.0_dp*step)
      call response_group_errors(finite_difference, analytic_action, current_error)
      minimum_error = min(minimum_error, current_error)
    end do

    einstein_error_n = maxval(abs(base_response%diffusion_n &
      - (BOLTZMANN_CONSTANT*base_state%temperature/ELEMENTARY_CHARGE) &
      * base_response%mobility_n))
    einstein_error_p = maxval(abs(base_response%diffusion_p &
      - (BOLTZMANN_CONSTANT*base_state%temperature/ELEMENTARY_CHARGE) &
      * base_response%mobility_p))

    passed = maxval(minimum_error) < 2.0e-7_dp .and. &
             neutrality_error < 1.0e-14_dp .and. mass_action_error < 1.0e-14_dp .and. &
             einstein_error_n < 1.0e-15_dp .and. einstein_error_p < 1.0e-15_dp .and. &
             base_response%n > 0.0_dp .and. base_response%p > 0.0_dp .and. &
             base_response%joule_heat >= 0.0_dp

    write(*,'(a)') 'verification_id=G02_material_point'
    call print_metric('stress_tangent_error', minimum_error(1))
    call print_metric('electric_displacement_tangent_error', minimum_error(2))
    call print_metric('charge_tangent_error', minimum_error(3))
    call print_metric('electron_current_tangent_error', minimum_error(4))
    call print_metric('hole_current_tangent_error', minimum_error(5))
    call print_metric('heat_flux_tangent_error', minimum_error(6))
    call print_metric('joule_heat_tangent_error', minimum_error(7))
    call print_metric('recombination_tangent_error', minimum_error(8))
    call print_metric('equilibrium_neutrality_error', neutrality_error)
    call print_metric('equilibrium_mass_action_error', mass_action_error)
    call print_metric('einstein_error_n', einstein_error_n)
    call print_metric('einstein_error_p', einstein_error_p)
    call print_metric('material_point_n', base_response%n)
    call print_metric('material_point_p', base_response%p)
    call print_metric('material_point_joule_heat', base_response%joule_heat)
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status, STATUS_VERIFICATION_FAILED, &
        'G02 material-point response or tangent check failed.')
    end if
  end subroutine verify_material_point_g02

  subroutine define_test_material(material)
    type(MaterialType), intent(out) :: material
    call define_tps_g08_reference_material(material)
    material%name='G02_ANISOTROPIC_TEST'
  end subroutine define_test_material

  subroutine define_test_state(state)
    type(MaterialPointStateType), intent(out) :: state
    state = MaterialPointStateType()
    state%strain = [1.0e-4_dp,-0.7e-4_dp,0.5e-4_dp,0.3e-4_dp,-0.2e-4_dp,0.4e-4_dp]
    state%electric_field = [1.2e5_dp,-0.8e5_dp,0.6e5_dp]
    state%temperature = 345.0_dp
    state%temperature_gradient = [1.0e3_dp,-2.0e3_dp,1.5e3_dp]
    state%eta_n = 0.12_dp
    state%eta_n_gradient = [1.1e3_dp,-0.7e3_dp,0.4e3_dp]
    state%eta_p = -0.18_dp
    state%eta_p_gradient = [-0.8e3_dp,0.9e3_dp,-0.5e3_dp]
  end subroutine define_test_state

  function input_scales() result(scales)
    real(dp) :: scales(INPUT_SIZE)
    scales(1:6) = 1.0e-4_dp
    scales(7:9) = 1.0e5_dp
    scales(10) = 10.0_dp
    scales(11:13) = 1.0e3_dp
    scales(14) = 0.1_dp
    scales(15:17) = 1.0e3_dp
    scales(18) = 0.1_dp
    scales(19:21) = 1.0e3_dp
  end function input_scales

  subroutine perturb_state(base, direction, factor, perturbed)
    type(MaterialPointStateType), intent(in) :: base
    real(dp), intent(in) :: direction(INPUT_SIZE), factor
    type(MaterialPointStateType), intent(out) :: perturbed

    perturbed = base
    perturbed%strain = base%strain + factor*direction(1:6)
    perturbed%electric_field = base%electric_field + factor*direction(7:9)
    perturbed%temperature = base%temperature + factor*direction(10)
    perturbed%temperature_gradient = base%temperature_gradient + factor*direction(11:13)
    perturbed%eta_n = base%eta_n + factor*direction(14)
    perturbed%eta_n_gradient = base%eta_n_gradient + factor*direction(15:17)
    perturbed%eta_p = base%eta_p + factor*direction(18)
    perturbed%eta_p_gradient = base%eta_p_gradient + factor*direction(19:21)
  end subroutine perturb_state

  subroutine response_group_errors(finite_difference, analytic_action, errors)
    real(dp), intent(in) :: finite_difference(RESPONSE_SIZE), analytic_action(RESPONSE_SIZE)
    real(dp), intent(out) :: errors(NUMBER_OF_GROUPS)
    errors(1) = relative_error(finite_difference(1:6),analytic_action(1:6))
    errors(2) = relative_error(finite_difference(7:9),analytic_action(7:9))
    errors(3) = relative_error(finite_difference(10:10),analytic_action(10:10))
    errors(4) = relative_error(finite_difference(11:13),analytic_action(11:13))
    errors(5) = relative_error(finite_difference(14:16),analytic_action(14:16))
    errors(6) = relative_error(finite_difference(17:19),analytic_action(17:19))
    errors(7) = relative_error(finite_difference(20:20),analytic_action(20:20))
    errors(8) = relative_error(finite_difference(21:21),analytic_action(21:21))
  end subroutine response_group_errors

  real(dp) function relative_error(value_a, value_b)
    real(dp), intent(in) :: value_a(:), value_b(:)
    real(dp) :: denominator
    denominator = max(1.0e-30_dp, sqrt(sum(value_a*value_a)), sqrt(sum(value_b*value_b)))
    relative_error = sqrt(sum((value_a-value_b)**2)) / denominator
  end function relative_error

  subroutine print_metric(name, value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name) // '=', value
  end subroutine print_metric

end module material_tangent_check_m

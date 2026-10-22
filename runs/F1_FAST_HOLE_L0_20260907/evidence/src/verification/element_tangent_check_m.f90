module element_tangent_check_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : NDOF_NODE, F_U1, F_U2, F_U3, F_PHI, F_ETA_N, F_ETA_P, &
    F_THETA, local_dof
  use mesh_types_m, only : HEX20_NNODE
  use hex20_shape_m, only : hex20_natural_coordinates
  use dof_map_m, only : HEX20_NDOF
  use model_types_m, only : MaterialType, LoadStateType
  use element_kernel_m, only : ElementRequestType, ElementDiagnosticsType, evaluate_hex20_element
  use material_tangent_check_m, only : define_test_material
  implicit none
  private

  integer, parameter :: NUMBER_OF_GROUPS = 5
  public :: verify_element_g03

contains

  subroutine verify_element_g03(status)
    type(SolverStatusType), intent(out) :: status
    type(MaterialType) :: material
    type(LoadStateType) :: load
    type(ElementRequestType) :: request
    type(ElementDiagnosticsType) :: diagnostics
    type(SolverStatusType) :: local_status
    real(dp) :: natural_nodes(3,HEX20_NNODE), coordinates(3,HEX20_NNODE)
    real(dp) :: state_e(HEX20_NDOF), direction(HEX20_NDOF), raw(HEX20_NDOF)
    real(dp) :: state_plus(HEX20_NDOF), state_minus(HEX20_NDOF)
    real(dp) :: residual(HEX20_NDOF), residual_plus(HEX20_NDOF), residual_minus(HEX20_NDOF)
    real(dp) :: jacobian(HEX20_NDOF,HEX20_NDOF), unused_jacobian(HEX20_NDOF,HEX20_NDOF)
    real(dp) :: analytic_action(HEX20_NDOF), finite_difference(HEX20_NDOF)
    real(dp) :: minimum_error(NUMBER_OF_GROUPS), current_error(NUMBER_OF_GROUPS)
    real(dp) :: step_values(8), step, electron_balance, hole_balance
    integer :: index, node, field, step_index
    integer :: indices_u(3*HEX20_NNODE), indices_phi(HEX20_NNODE)
    integer :: indices_n(HEX20_NNODE), indices_p(HEX20_NNODE), indices_T(HEX20_NNODE)
    logical :: passed

    call clear_status(status)
    call define_test_material(material)
    material%reference_p = 2.0e20_dp
    call hex20_natural_coordinates(natural_nodes)
    coordinates = 0.5_dp * (natural_nodes + 1.0_dp)
    call define_element_state(coordinates, state_e)
    load%voltage_factor = 1.0_dp
    load%joule_factor = 0.81_dp
    request%need_jacobian = .true.
    call evaluate_hex20_element(material, coordinates, state_e, load, request, &
      residual, jacobian, diagnostics, local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
      return
    end if

    call build_field_indices(indices_u, indices_phi, indices_n, indices_p, indices_T)
    do index = 1, HEX20_NDOF
      raw(index) = sin(0.23_dp*real(index,dp)) + cos(0.41_dp*real(index,dp))
    end do
    raw = raw / sqrt(sum(raw*raw))
    do node = 1, HEX20_NNODE
      do field = 1, NDOF_NODE
        index = local_dof(node,field)
        select case (field)
        case (F_U1,F_U2,F_U3)
          direction(index) = raw(index)*1.0e-5_dp
        case (F_PHI)
          direction(index) = raw(index)*1.0e3_dp
        case (F_ETA_N,F_ETA_P)
          direction(index) = raw(index)*1.0e-2_dp
        case (F_THETA)
          direction(index) = raw(index)*1.0_dp
        end select
      end do
    end do
    analytic_action = matmul(jacobian,direction)

    step_values = [1.0e-2_dp,3.0e-3_dp,1.0e-3_dp,3.0e-4_dp,1.0e-4_dp, &
                   3.0e-5_dp,1.0e-5_dp,3.0e-6_dp]
    minimum_error = huge(1.0_dp)
    request%need_jacobian = .false.
    do step_index = 1, size(step_values)
      step = step_values(step_index)
      state_plus = state_e + step*direction
      state_minus = state_e - step*direction
      call evaluate_hex20_element(material, coordinates, state_plus, load, request, &
        residual_plus, unused_jacobian, diagnostics, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
        return
      end if
      call evaluate_hex20_element(material, coordinates, state_minus, load, request, &
        residual_minus, unused_jacobian, diagnostics, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
        return
      end if
      finite_difference = (residual_plus-residual_minus)/(2.0_dp*step)
      current_error(1) = relative_error(finite_difference(indices_u),analytic_action(indices_u))
      current_error(2) = relative_error(finite_difference(indices_phi),analytic_action(indices_phi))
      current_error(3) = relative_error(finite_difference(indices_n),analytic_action(indices_n))
      current_error(4) = relative_error(finite_difference(indices_p),analytic_action(indices_p))
      current_error(5) = relative_error(finite_difference(indices_T),analytic_action(indices_T))
      minimum_error = min(minimum_error,current_error)
    end do

    electron_balance = abs(sum(residual(indices_n))) / max(1.0_dp,vector_norm(residual(indices_n)))
    hole_balance = abs(sum(residual(indices_p))) / max(1.0_dp,vector_norm(residual(indices_p)))
    passed = maxval(minimum_error) < 5.0e-7_dp .and. &
      electron_balance < 5.0e-12_dp .and. hole_balance < 5.0e-12_dp .and. &
      abs(diagnostics%volume-1.0_dp) < 5.0e-12_dp .and. &
      diagnostics%min_n > 0.0_dp .and. diagnostics%min_p > 0.0_dp .and. &
      diagnostics%min_temperature > 0.0_dp

    write(*,'(a)') 'verification_id=G03_element_kernel'
    call print_metric('mechanical_row_tangent_error',minimum_error(1))
    call print_metric('potential_row_tangent_error',minimum_error(2))
    call print_metric('electron_row_tangent_error',minimum_error(3))
    call print_metric('hole_row_tangent_error',minimum_error(4))
    call print_metric('thermal_row_tangent_error',minimum_error(5))
    call print_metric('electron_constant_test_balance',electron_balance)
    call print_metric('hole_constant_test_balance',hole_balance)
    call print_metric('element_volume',diagnostics%volume)
    call print_metric('element_min_n',diagnostics%min_n)
    call print_metric('element_min_p',diagnostics%min_p)
    call print_metric('element_min_temperature',diagnostics%min_temperature)
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status, STATUS_VERIFICATION_FAILED, &
        'G03 element residual, Jacobian, or conservation check failed.')
    end if
  end subroutine verify_element_g03

  subroutine define_element_state(coordinates, state_e)
    real(dp), intent(in) :: coordinates(3,HEX20_NNODE)
    real(dp), intent(out) :: state_e(HEX20_NDOF)
    real(dp) :: x, y, z
    integer :: node

    state_e = 0.0_dp
    do node = 1, HEX20_NNODE
      x = coordinates(1,node)
      y = coordinates(2,node)
      z = coordinates(3,node)
      state_e(local_dof(node,F_U1)) = 1.0e-4_dp*x + 0.2e-4_dp*y - 0.1e-4_dp*z
      state_e(local_dof(node,F_U2)) = -0.3e-4_dp*x + 0.8e-4_dp*y + 0.2e-4_dp*z
      state_e(local_dof(node,F_U3)) = 0.1e-4_dp*x - 0.2e-4_dp*y + 0.6e-4_dp*z
      state_e(local_dof(node,F_PHI)) = -1.2e5_dp*x + 0.8e5_dp*y - 0.6e5_dp*z
      state_e(local_dof(node,F_ETA_N)) = 0.12_dp + 0.05_dp*x - 0.03_dp*y + 0.02_dp*z
      state_e(local_dof(node,F_ETA_P)) = -0.18_dp - 0.04_dp*x + 0.06_dp*y - 0.01_dp*z
      state_e(local_dof(node,F_THETA)) = 45.0_dp + 3.0_dp*x - 2.0_dp*y + 1.5_dp*z
    end do
  end subroutine define_element_state

  subroutine build_field_indices(indices_u, indices_phi, indices_n, indices_p, indices_T)
    integer, intent(out) :: indices_u(3*HEX20_NNODE), indices_phi(HEX20_NNODE)
    integer, intent(out) :: indices_n(HEX20_NNODE), indices_p(HEX20_NNODE)
    integer, intent(out) :: indices_T(HEX20_NNODE)
    integer :: node, component, counter

    counter = 0
    do node = 1, HEX20_NNODE
      do component = F_U1, F_U3
        counter = counter + 1
        indices_u(counter) = local_dof(node,component)
      end do
      indices_phi(node) = local_dof(node,F_PHI)
      indices_n(node) = local_dof(node,F_ETA_N)
      indices_p(node) = local_dof(node,F_ETA_P)
      indices_T(node) = local_dof(node,F_THETA)
    end do
  end subroutine build_field_indices

  real(dp) function relative_error(value_a, value_b)
    real(dp), intent(in) :: value_a(:), value_b(:)
    real(dp) :: denominator
    denominator = max(1.0e-30_dp,vector_norm(value_a),vector_norm(value_b))
    relative_error = vector_norm(value_a-value_b)/denominator
  end function relative_error

  real(dp) function vector_norm(vector)
    real(dp), intent(in) :: vector(:)
    vector_norm = sqrt(sum(vector*vector))
  end function vector_norm

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

end module element_tangent_check_m

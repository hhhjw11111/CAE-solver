module element_kernel_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_GEOMETRY_ERROR, STATUS_INVALID_ARGUMENT
  use field_layout_m, only : F_U1, F_U2, F_U3, F_PHI, F_ETA_N, F_ETA_P, F_THETA, local_dof
  use mesh_types_m, only : HEX20_NNODE
  use hex20_quadrature_m, only : HEX20_NGAUSS, gauss_rule_3x3x3
  use hex20_shape_m, only : hex20_geometry
  use transport_supg_m, only : hex20_hessians,add_transport_supg
  use dof_map_m, only : HEX20_NDOF
  use model_types_m, only : MaterialType, LoadStateType
  use material_point_m, only : MATERIAL_INPUT_SIZE, MATERIAL_RESPONSE_SIZE, &
    MaterialPointStateType, MaterialPointResponseType, MaterialPointTangentType, &
    evaluate_material_point, pack_material_response, material_tangent_matrix
  implicit none
  private

  type, public :: ElementRequestType
    logical :: need_jacobian = .true.
  end type ElementRequestType

  type, public :: ElementDiagnosticsType
    real(dp) :: volume = 0.0_dp
    real(dp) :: min_temperature = huge(1.0_dp)
    real(dp) :: max_temperature = -huge(1.0_dp)
    real(dp) :: min_n = huge(1.0_dp)
    real(dp) :: max_n = -huge(1.0_dp)
    real(dp) :: min_p = huge(1.0_dp)
    real(dp) :: max_p = -huge(1.0_dp)
    real(dp) :: max_joule_heat = 0.0_dp
  end type ElementDiagnosticsType

  public :: evaluate_hex20_element, build_operators, unpack_point_state

contains

  subroutine evaluate_hex20_element(material, coordinates, state_e, load, request, &
                                    residual_e, jacobian_e, diagnostics, status)
    type(MaterialType), intent(in) :: material
    real(dp), intent(in) :: coordinates(3,HEX20_NNODE)
    real(dp), intent(in) :: state_e(HEX20_NDOF)
    type(LoadStateType), intent(in) :: load
    type(ElementRequestType), intent(in) :: request
    real(dp), intent(out) :: residual_e(HEX20_NDOF)
    real(dp), intent(out) :: jacobian_e(HEX20_NDOF,HEX20_NDOF)
    type(ElementDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    real(dp) :: points(3,HEX20_NGAUSS), weights(HEX20_NGAUSS)
    real(dp) :: shape(HEX20_NNODE), derivative_global(3,HEX20_NNODE)
    real(dp) :: determinant, integration_weight
    real(dp) :: hess(3,3,HEX20_NNODE),metric(3,3)
    real(dp) :: state_operator(MATERIAL_INPUT_SIZE,HEX20_NDOF)
    real(dp) :: test_operator(HEX20_NDOF,MATERIAL_RESPONSE_SIZE)
    real(dp) :: independent_state(MATERIAL_INPUT_SIZE)
    real(dp) :: response_vector(MATERIAL_RESPONSE_SIZE)
    real(dp) :: tangent_matrix(MATERIAL_RESPONSE_SIZE,MATERIAL_INPUT_SIZE)
    real(dp) :: log_weights(2)
    type(MaterialPointStateType) :: point_state
    type(MaterialPointResponseType) :: point_response
    type(MaterialPointTangentType) :: point_tangent
    type(SolverStatusType) :: local_status
    integer :: point

    call clear_status(status)
    residual_e = 0.0_dp
    jacobian_e = 0.0_dp
    diagnostics = ElementDiagnosticsType()
    if(material%transport_log_weighted.and.(material%transport_supg.or.material%linear_carriers)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Log-weighted transport excludes SUPG and linear carriers.')
      return
    end if
    call gauss_rule_3x3x3(points, weights)

    do point = 1, HEX20_NGAUSS
      call hex20_geometry(coordinates, points(1,point), points(2,point), points(3,point), &
        shape, derivative_global, determinant, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, STATUS_GEOMETRY_ERROR, &
          'Element geometry failed at Gauss point ' // integer_to_string(point) // '.')
        return
      end if
      integration_weight = weights(point) * determinant
      call build_operators(shape, derivative_global, state_operator, test_operator)
      independent_state = matmul(state_operator, state_e)
      if(material%transport_log_weighted) then
        if(any(.not.ieee_is_finite(independent_state([14,18]))).or. &
            any(abs(independent_state([14,18]))>=log(huge(1.0_dp)))) then
          call set_status(status,STATUS_INVALID_ARGUMENT,'Log-weighted transport exponent is out of range.')
          return
        end if
        call log_transport_tests(shape,derivative_global,independent_state,test_operator,log_weights)
      end if
      call unpack_point_state(material, independent_state, point_state)
      call evaluate_material_point(material, point_state, load, point_response, &
        point_tangent, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, local_status%code, &
          'Material point failed at Gauss point ' // integer_to_string(point) // &
          ': ' // trim(local_status%message))
        return
      end if
      call pack_material_response(point_response, response_vector)
      residual_e = residual_e + integration_weight * matmul(test_operator, response_vector)
      if (request%need_jacobian) then
        call material_tangent_matrix(point_tangent, tangent_matrix)
        jacobian_e = jacobian_e + integration_weight * &
          matmul(test_operator, matmul(tangent_matrix, state_operator))
        if(material%transport_log_weighted) then
          call add_log_test_derivative(shape,derivative_global,test_operator,response_vector, &
            log_weights,integration_weight,jacobian_e)
        end if
      end if
      if (material%transport_supg) then
        call hex20_hessians(coordinates,points(:,point),hess,metric,local_status)
        if (.not.local_status%is_ok()) then
          status=local_status;return
        end if
        call add_transport_supg(material,shape,derivative_global,hess,metric,state_e, &
          integration_weight,request%need_jacobian,residual_e,jacobian_e,local_status)
        if (.not.local_status%is_ok()) then
          status=local_status;return
        end if
      end if
      diagnostics%volume = diagnostics%volume + integration_weight
      diagnostics%min_temperature = min(diagnostics%min_temperature, point_state%temperature)
      diagnostics%max_temperature = max(diagnostics%max_temperature, point_state%temperature)
      diagnostics%min_n = min(diagnostics%min_n, point_response%n)
      diagnostics%max_n = max(diagnostics%max_n, point_response%n)
      diagnostics%min_p = min(diagnostics%min_p, point_response%p)
      diagnostics%max_p = max(diagnostics%max_p, point_response%p)
      diagnostics%max_joule_heat = max(diagnostics%max_joule_heat, point_response%joule_heat)
    end do
  end subroutine evaluate_hex20_element

  subroutine log_transport_tests(shape,dn,q,test,weights)
    real(dp), intent(in) :: shape(20),dn(3,20),q(21)
    real(dp), intent(inout) :: test(140,21)
    real(dp), intent(out) :: weights(2)
    integer :: i,nrow,prow
    weights=exp(-q([14,18]))
    do i=1,20
      nrow=local_dof(i,F_ETA_N);prow=local_dof(i,F_ETA_P)
      test(nrow,11:13)=weights(1)*(dn(:,i)-shape(i)*q(15:17))
      test(prow,14:16)=weights(2)*(dn(:,i)-shape(i)*q(19:21))
      test(nrow,21)=weights(1)*ELEMENTARY_CHARGE*shape(i)
      test(prow,21)=-weights(2)*ELEMENTARY_CHARGE*shape(i)
    end do
  end subroutine log_transport_tests

  subroutine add_log_test_derivative(shape,dn,test,response,weights,integration_weight,k)
    real(dp), intent(in) :: shape(20),dn(3,20),test(140,21),response(21),weights(2),integration_weight
    real(dp), intent(inout) :: k(140,140)
    real(dp) :: base_n,base_p
    integer :: i,j,nrow,prow
    do i=1,20
      nrow=local_dof(i,F_ETA_N);prow=local_dof(i,F_ETA_P)
      base_n=dot_product(test(nrow,:),response);base_p=dot_product(test(prow,:),response)
      do j=1,20
        k(nrow,local_dof(j,F_ETA_N))=k(nrow,local_dof(j,F_ETA_N))-integration_weight*( &
          shape(j)*base_n+weights(1)*shape(i)*dot_product(dn(:,j),response(11:13)))
        k(prow,local_dof(j,F_ETA_P))=k(prow,local_dof(j,F_ETA_P))-integration_weight*( &
          shape(j)*base_p+weights(2)*shape(i)*dot_product(dn(:,j),response(14:16)))
      end do
    end do
  end subroutine add_log_test_derivative

  subroutine build_operators(shape, derivative_global, state_operator, test_operator)
    real(dp), intent(in) :: shape(HEX20_NNODE), derivative_global(3,HEX20_NNODE)
    real(dp), intent(out) :: state_operator(MATERIAL_INPUT_SIZE,HEX20_NDOF)
    real(dp), intent(out) :: test_operator(HEX20_NDOF,MATERIAL_RESPONSE_SIZE)
    integer :: node, column_u1, column_u2, column_u3, column_phi
    integer :: column_eta_n, column_eta_p, column_theta
    real(dp) :: dx, dy, dz

    state_operator = 0.0_dp
    test_operator = 0.0_dp
    do node = 1, HEX20_NNODE
      dx = derivative_global(1,node)
      dy = derivative_global(2,node)
      dz = derivative_global(3,node)
      column_u1 = local_dof(node,F_U1)
      column_u2 = local_dof(node,F_U2)
      column_u3 = local_dof(node,F_U3)
      column_phi = local_dof(node,F_PHI)
      column_eta_n = local_dof(node,F_ETA_N)
      column_eta_p = local_dof(node,F_ETA_P)
      column_theta = local_dof(node,F_THETA)

      state_operator(1,column_u1) = dx
      state_operator(5,column_u1) = dz
      state_operator(6,column_u1) = dy
      state_operator(2,column_u2) = dy
      state_operator(4,column_u2) = dz
      state_operator(6,column_u2) = dx
      state_operator(3,column_u3) = dz
      state_operator(4,column_u3) = dy
      state_operator(5,column_u3) = dx
      state_operator(7:9,column_phi) = -derivative_global(:,node)
      state_operator(14,column_eta_n) = shape(node)
      state_operator(15:17,column_eta_n) = derivative_global(:,node)
      state_operator(18,column_eta_p) = shape(node)
      state_operator(19:21,column_eta_p) = derivative_global(:,node)
      state_operator(10,column_theta) = shape(node)
      state_operator(11:13,column_theta) = derivative_global(:,node)

      test_operator(column_u1,1) = dx
      test_operator(column_u1,5) = dz
      test_operator(column_u1,6) = dy
      test_operator(column_u2,2) = dy
      test_operator(column_u2,4) = dz
      test_operator(column_u2,6) = dx
      test_operator(column_u3,3) = dz
      test_operator(column_u3,4) = dy
      test_operator(column_u3,5) = dx
      test_operator(column_phi,7:9) = derivative_global(:,node)
      test_operator(column_phi,10) = shape(node)
      test_operator(column_eta_n,11:13) = derivative_global(:,node)
      test_operator(column_eta_n,21) = ELEMENTARY_CHARGE * shape(node)
      test_operator(column_eta_p,14:16) = derivative_global(:,node)
      test_operator(column_eta_p,21) = -ELEMENTARY_CHARGE * shape(node)
      test_operator(column_theta,17:19) = -derivative_global(:,node)
      test_operator(column_theta,20) = -shape(node)
    end do
  end subroutine build_operators

  subroutine unpack_point_state(material, independent_state, point_state)
    type(MaterialType), intent(in) :: material
    real(dp), intent(in) :: independent_state(MATERIAL_INPUT_SIZE)
    type(MaterialPointStateType), intent(out) :: point_state

    point_state = MaterialPointStateType()
    point_state%strain = independent_state(1:6)
    point_state%electric_field = independent_state(7:9)
    point_state%temperature = material%reference_temperature + independent_state(10)
    point_state%temperature_gradient = independent_state(11:13)
    point_state%eta_n = independent_state(14)
    point_state%eta_n_gradient = independent_state(15:17)
    point_state%eta_p = independent_state(18)
    point_state%eta_p_gradient = independent_state(19:21)
  end subroutine unpack_point_state

  function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

end module element_kernel_m

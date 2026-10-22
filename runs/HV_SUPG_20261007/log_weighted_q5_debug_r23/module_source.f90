! Verification-only prototype. The production kernel and physical material
! response are linked from an immutable build; no production mode is changed.
module log_weighted_probe_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use status_m, only : SolverStatusType,set_status,STATUS_INVALID_ARGUMENT
  use model_types_m, only : MaterialType,LoadStateType
  use field_layout_m, only : F_ETA_N,F_ETA_P,local_dof
  use hex20_shape_m, only : hex20_geometry
  use hex20_quadrature_m, only : gauss_rule_3x3x3
  use element_kernel_m, only : ElementRequestType,ElementDiagnosticsType, &
    evaluate_hex20_element,build_operators,unpack_point_state
  use material_point_m, only : MaterialPointStateType,MaterialPointResponseType, &
    MaterialPointTangentType,evaluate_material_point,pack_material_response,material_tangent_matrix
  implicit none
  private
  public :: evaluate_log_hex20
contains
  subroutine evaluate_log_hex20(material,coordinates,state_e,load,request,r,k,diagnostics,status)
    type(MaterialType), intent(in) :: material
    real(dp), intent(in) :: coordinates(3,20),state_e(140)
    type(LoadStateType), intent(in) :: load
    type(ElementRequestType), intent(in) :: request
    real(dp), intent(out) :: r(140),k(140,140)
    type(ElementDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    type(MaterialPointStateType) :: state
    type(MaterialPointResponseType) :: response
    type(MaterialPointTangentType) :: tangent
    real(dp) :: points(3,27),weights(27),n(20),dn(3,20),det,weight
    real(dp) :: sop(21,140),top(140,21),q(21),v(21),a(21,21),eta,g(3),b(3),test(21)
    real(dp) :: factor,reaction_sign,base
    integer :: gp,carrier,field,eta_index,current_index,i,j,row,col
    if(material%transport_supg.or.material%linear_carriers) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Log-weighted probe excludes SUPG and LQ.');return
    end if
    ! Original mechanics, electrostatics and heat rows are retained bit for bit.
    call evaluate_hex20_element(material,coordinates,state_e,load,request,r,k,diagnostics,status)
    if(.not.status%is_ok()) return
    r(F_ETA_N::7)=0.0_dp;r(F_ETA_P::7)=0.0_dp
    k(F_ETA_N::7,:)=0.0_dp;k(F_ETA_P::7,:)=0.0_dp
    call gauss_rule_3x3x3(points,weights)
    do gp=1,27
      call hex20_geometry(coordinates,points(1,gp),points(2,gp),points(3,gp),n,dn,det,status)
      if(.not.status%is_ok()) return
      weight=det*weights(gp)
      call build_operators(n,dn,sop,top)
      q=matmul(sop,state_e)
      call unpack_point_state(material,q,state)
      call evaluate_material_point(material,state,load,response,tangent,status)
      if(.not.status%is_ok()) return
      call pack_material_response(response,v)
      if(request%need_jacobian) call material_tangent_matrix(tangent,a)
      do carrier=1,2
        field=F_ETA_N;eta_index=14;current_index=11;reaction_sign=1.0_dp
        if(carrier==2) then
          field=F_ETA_P;eta_index=18;current_index=14;reaction_sign=-1.0_dp
        end if
        eta=q(eta_index);g=q(eta_index+1:eta_index+3)
        factor=exp(-eta)
        do i=1,20
          row=local_dof(i,field);b=dn(:,i)-n(i)*g
          test=0.0_dp
          test(current_index:current_index+2)=factor*b
          test(21)=factor*reaction_sign*ELEMENTARY_CHARGE*n(i)
          base=dot_product(b,v(current_index:current_index+2))+ &
            reaction_sign*ELEMENTARY_CHARGE*n(i)*v(21)
          r(row)=r(row)+weight*factor*base
          if(request%need_jacobian) then
            k(row,:)=k(row,:)+weight*matmul(matmul(test,a),sop)
            ! State-dependent test derivative: d[e^-eta(grad N-N grad eta)].
            do j=1,20
              col=local_dof(j,field)
              k(row,col)=k(row,col)-weight*factor*(n(j)*base+ &
                n(i)*dot_product(dn(:,j),v(current_index:current_index+2)))
            end do
          end if
        end do
      end do
    end do
  end subroutine evaluate_log_hex20
end module log_weighted_probe_m

! Consistent residual-based SUPG for the existing log-density drift-diffusion PDE.
! This is NOT a Scharfetter-Gummel flux or an energy-model modification.
module transport_supg_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use field_layout_m, only : F_PHI,F_ETA_N,F_ETA_P,F_THETA,local_dof
  use model_types_m, only : MaterialType
  use mesh_types_m, only : HEX20_NNODE
  use dof_map_m, only : HEX20_NDOF
  use hex20_shape_m, only : hex20_shape,hex20_geometry,inverse_3x3
  use transport_m, only : evaluate_transport
  use status_m, only : SolverStatusType,set_status,STATUS_INVALID_ARGUMENT
  implicit none
  private
  public :: hex20_hessians, add_transport_supg
contains
  subroutine hex20_hessians(x,xi,hess,metric,status)
    real(dp), intent(in) :: x(3,20),xi(3)
    real(dp), intent(out) :: hess(3,3,20),metric(3,3)
    type(SolverStatusType), intent(out) :: status
    real(dp) :: n(20),dn(3,20),dg(3,20),plus(3,20),minus(3,20)
    real(dp) :: hn(3,3,20),hx(3,3,3),j(3,3),ji(3,3),det,r(3),tmp(3,3)
    integer :: a,b,k
    call hex20_geometry(x,xi(1),xi(2),xi(3),n,dg,det,status)
    if (.not.status%is_ok()) return
    call hex20_shape(xi(1),xi(2),xi(3),n,dn)
    j=matmul(dn,transpose(x))
    call inverse_3x3(j,ji,det)
    metric=matmul(ji,transpose(ji))
    ! Each first derivative is at most quadratic in each natural coordinate.
    ! Centered difference with unit offset is its EXACT polynomial derivative.
    do b=1,3
      r=xi;r(b)=r(b)+1.0_dp
      call hex20_shape(r(1),r(2),r(3),n,plus)
      r=xi;r(b)=r(b)-1.0_dp
      call hex20_shape(r(1),r(2),r(3),n,minus)
      hn(:,b,:)=(plus-minus)*0.5_dp
    end do
    hx=0.0_dp
    do k=1,3
      do a=1,20
        hx(:,:,k)=hx(:,:,k)+x(k,a)*hn(:,:,a)
      end do
    end do
    do a=1,20
      tmp=hn(:,:,a)
      do k=1,3
        tmp=tmp-dg(k,a)*hx(:,:,k)
      end do
      hess(:,:,a)=matmul(ji,matmul(tmp,transpose(ji)))
    end do
  end subroutine hex20_hessians

  subroutine add_transport_supg(m,n,dn,hess,metric,x,w,need_jac,r,k,status)
    type(MaterialType), intent(in) :: m
    real(dp), intent(in) :: n(20),dn(3,20),hess(3,3,20),metric(3,3)
    real(dp), intent(in) :: x(140),w
    logical, intent(in) :: need_jac
    real(dp), intent(inout) :: r(140),k(140,140)
    type(SolverStatusType), intent(out) :: status
    real(dp) :: mu(3,3,2),d(3,3,2),mut(3,3,2),dt(3,3,2),mutt(3,3),dtt(3,3)
    real(dp) :: e(3),gt(3),g(3),hp(3,3),hg(3,3),t,eta,rho,beta,sgn,ref
    real(dp) :: v(3),vt(3),divv,strong,ds,drho,tau,den,diff,trg,adv(20)
    real(dp) :: de(3),dgt(3),dg(3),dhp(3,3),dhg(3,3),dtemp,deta,dv(3),dvt(3)
    real(dp) :: dm(3,3),dd(3,3),dmt(3,3),ddt(3,3),ddivv,ddiff,dtau,dsadv(20)
    real(dp) :: velocity2,advmetric,epsilon,gdir,dgdir,depsilon,dvelocity2,dadvmetric
    integer :: c,a,b,f,col,rows(20),field
    if(m%linear_carriers .or. m%supg_scale<=0.0_dp) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'SUPG requires nonlinear log densities and a positive scale.')
      return
    end if
    e=0.0_dp;gt=0.0_dp;hp=0.0_dp;t=m%reference_temperature
    do a=1,20
      e=e-dn(:,a)*x(local_dof(a,F_PHI))
      hp=hp+hess(:,:,a)*x(local_dof(a,F_PHI))
      t=t+n(a)*x(local_dof(a,F_THETA))
      gt=gt+dn(:,a)*x(local_dof(a,F_THETA))
    end do
    call evaluate_transport(m,t,mu(:,:,1),mu(:,:,2),d(:,:,1),d(:,:,2), &
      mut(:,:,1),mut(:,:,2),dt(:,:,1),dt(:,:,2),status)
    if (.not.status%is_ok()) return
    trg=metric(1,1)+metric(2,2)+metric(3,3)
    do c=1,2
      field=F_ETA_N;sgn=1.0_dp;ref=m%reference_n;beta=m%beta_n
      if(c==2) then
        field=F_ETA_P;sgn=-1.0_dp;ref=m%reference_p;beta=m%beta_p
      end if
      eta=0.0_dp;g=0.0_dp;hg=0.0_dp
      do a=1,20
        rows(a)=local_dof(a,field)
        eta=eta+n(a)*x(rows(a));g=g+dn(:,a)*x(rows(a));hg=hg+hess(:,:,a)*x(rows(a))
      end do
      rho=ref*exp(eta)
      v=matmul(mu(:,:,c),e);vt=matmul(mut(:,:,c),e)
      divv=-sum(mu(:,:,c)*hp)+dot_product(vt,gt)
      strong=dot_product(v,g)+divv+sgn*(sum(d(:,:,c)*hg)+ &
        dot_product(g,matmul(d(:,:,c),g))+dot_product(gt,matmul(dt(:,:,c),g)))
      diff=(d(1,1,c)+d(2,2,c)+d(3,3,c))/3.0_dp
      if(diff<=0.0_dp .or.rho<=0.0_dp) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'SUPG requires positive diffusion and density.');return
      end if
      ! Smooth metric tau, quadratic-element diffusion constant 4.
      velocity2=dot_product(v,v);advmetric=dot_product(v,matmul(metric,v))
      epsilon=1.0e-12_dp*diff*diff*trg
      gdir=trg
      if(m%supg_directional_tau) gdir=(advmetric+epsilon*trg/3.0_dp)/(velocity2+epsilon)
      ! Directional metric avoids suppressing streamwise stabilization solely
      ! because a transverse boundary layer has very small thickness elements.
      ! epsilon gives a smooth, differentiable zero-velocity limit.
      den=advmetric+(4.0_dp*diff*gdir)**2
      if(den<=tiny(1.0_dp)) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'SUPG requires nonzero transport coefficients.');return
      end if
      tau=m%supg_scale/sqrt(den)
      adv=matmul(transpose(dn),v)
      r(rows)=r(rows)+w*sgn*ELEMENTARY_CHARGE*rho*tau*adv*strong
      if(.not.need_jac) cycle
      mutt=beta*(beta+1.0_dp)*mu(:,:,c)/(t*t)
      dtt=beta*(beta-1.0_dp)*d(:,:,c)/(t*t)
      if(m%frozen_transport) then
        mutt=0.0_dp;dtt=0.0_dp
      end if
      do b=1,20
        do f=1,3
          de=0.0_dp;dgt=0.0_dp;dg=0.0_dp;dhp=0.0_dp;dhg=0.0_dp;dtemp=0.0_dp;deta=0.0_dp
          select case(f)
          case(1)
            col=local_dof(b,F_PHI);de=-dn(:,b);dhp=hess(:,:,b)
          case(2)
            col=rows(b);deta=n(b);dg=dn(:,b);dhg=hess(:,:,b)
          case(3)
            col=local_dof(b,F_THETA);dtemp=n(b);dgt=dn(:,b)
          end select
          dm=mut(:,:,c)*dtemp;dd=dt(:,:,c)*dtemp
          dmt=mutt*dtemp;ddt=dtt*dtemp
          dv=matmul(dm,e)+matmul(mu(:,:,c),de)
          dvt=matmul(dmt,e)+matmul(mut(:,:,c),de)
          ddivv=-sum(dm*hp)-sum(mu(:,:,c)*dhp)+dot_product(dvt,gt)+dot_product(vt,dgt)
          ds=dot_product(dv,g)+dot_product(v,dg)+ddivv+sgn*( &
            sum(dd*hg)+sum(d(:,:,c)*dhg)+dot_product(dg,matmul(d(:,:,c),g))+ &
            dot_product(g,matmul(dd,g))+dot_product(g,matmul(d(:,:,c),dg))+ &
            dot_product(dgt,matmul(dt(:,:,c),g))+dot_product(gt,matmul(ddt,g))+ &
            dot_product(gt,matmul(dt(:,:,c),dg)))
          ddiff=(dd(1,1)+dd(2,2)+dd(3,3))/3.0_dp
          dgdir=0.0_dp
          if(m%supg_directional_tau) then
            depsilon=2.0_dp*epsilon*ddiff/diff
            dvelocity2=2.0_dp*dot_product(v,dv)
            dadvmetric=2.0_dp*dot_product(dv,matmul(metric,v))
            dgdir=(dadvmetric+depsilon*trg/3.0_dp-gdir*(dvelocity2+depsilon))/(velocity2+epsilon)
          end if
          dtau=-tau*(dot_product(dv,matmul(metric,v))+16.0_dp*diff*ddiff*gdir*gdir+ &
            16.0_dp*diff*diff*gdir*dgdir)/den
          drho=rho*deta
          dsadv=matmul(transpose(dn),dv)
          k(rows,col)=k(rows,col)+w*sgn*ELEMENTARY_CHARGE*( &
            (drho*tau+rho*dtau)*adv*strong+rho*tau*(dsadv*strong+adv*ds))
        end do
      end do
    end do
  end subroutine add_transport_supg
end module transport_supg_m

module supg_trace_iso_check_m
  use kinds_m, only : dp
  use constants_m, only : BOLTZMANN_CONSTANT,ELEMENTARY_CHARGE
  use field_layout_m, only : F_PHI,F_ETA_N,F_ETA_P,F_THETA,local_dof
  use model_types_m, only : MaterialType,LoadStateType
  use hex20_shape_m, only : hex20_natural_coordinates
  use transport_supg_m, only : hex20_hessians
  use element_kernel_m, only : evaluate_hex20_element,ElementRequestType,ElementDiagnosticsType
  use material_tangent_check_m, only : define_test_material
  use element_tangent_check_m, only : verify_element_g03
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  implicit none
  private
  public :: verify_supg_trace_iso
contains
  subroutine verify_supg_trace_iso(status)
    type(SolverStatusType), intent(out) :: status
    real(dp) :: nat(3,20),x(3,20),h(3,3,20),metric(3,3),e,quad(3,3),v(20)
    real(dp) :: pe,err,last,cur,low_difference
    integer :: a,b,l,ip,c
    call clear_status(status)
    call verify_element_g03(status,.true.,.false.)
    if(.not.status%is_ok()) return
    call hex20_natural_coordinates(nat)
    x=0.5_dp*(nat+1.0_dp)
    do a=1,20
      x(1,a)=x(1,a)+0.13_dp*nat(2,a)*nat(3,a)
      x(2,a)=x(2,a)+0.08_dp*nat(1,a)*nat(3,a)
    end do
    call hex20_hessians(x,[0.21_dp,-0.37_dp,0.14_dp],h,metric,status)
    if(.not.status%is_ok()) return
    e=maxval(abs(sum(h,dim=3)))
    do b=1,3
      quad=0.0_dp
      do a=1,20
        quad=quad+x(b,a)*h(:,:,a)
      end do
      e=max(e,maxval(abs(quad)))
    end do
    write(*,'(a,es24.16e3)') 'supg_curved_geometry_linear_hessian_error=',e
    if(e>1.0e-11_dp) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG curved geometry Hessian failed.');return
    end if
    ! Independently known exponential concentration, nonzero constant current.
    ! Prescribed isothermal field, both carrier signs, production HEX20 kernel.
    do ip=1,6
      select case(ip)
      case(1);pe=0.1_dp
      case(2);pe=1.0_dp
      case(3);pe=10.0_dp
      case(4);pe=50.0_dp
      case(5);pe=200.0_dp
      case(6);pe=1000.0_dp
      end select
      do c=1,2
        last=huge(1.0_dp)
        do l=1,merge(4,3,ip==6)
          call carrier_bar(pe,8*2**l,c,.true.,err,cur,v,status)
          if(.not.status%is_ok()) return
          write(*,'(a,f6.1,a,i0,a,i0,a,es12.4,a,es12.4)') &
            'supg_bar_Pe=',pe,',carrier=',c,',elements=',8*2**l,',density_error=',err,',current_error=',cur
          if(err>last*1.05_dp) then
            call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG bar refinement failed.');return
          end if
          last=err
        end do
        if(err>1.0e-3_dp .or.cur>1.0e-2_dp) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG bar accuracy failed.');return
        end if
      end do
    end do
    call carrier_bar(1.0e-4_dp,16,1,.true.,err,cur,v,status)
    if(.not.status%is_ok()) return
    block
      real(dp) :: vg(20)
      call carrier_bar(1.0e-4_dp,16,1,.false.,err,cur,vg,status)
      if(.not.status%is_ok()) return
      low_difference=maxval(abs(v-vg))
    end block
    write(*,'(a,es24.16e3)') 'supg_low_pe_galerkin_difference=',low_difference
    if(low_difference>1.0e-7_dp) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG low Pe regression failed.');return
    end if
    write(*,'(a)') 'verification_id=TRACE_SUPG_ISO_HIGH_PE'
    write(*,'(a)') 'verification_result=PASS'
  end subroutine verify_supg_trace_iso

  subroutine carrier_bar(pe,ne,carrier,enabled,error,current_error,samples,status)
    real(dp), intent(in) :: pe
    integer, intent(in) :: ne,carrier
    logical, intent(in) :: enabled
    real(dp), intent(out) :: error,current_error,samples(20)
    type(SolverStatusType), intent(out) :: status
    type(MaterialType) :: m
    type(LoadStateType) :: load
    type(ElementRequestType) :: req
    type(ElementDiagnosticsType) :: diag
    real(dp) :: nat(3,20),coords(3,20),xe(140),re(140),ke(140,140)
    real(dp) :: eta(2*ne+1),res(2*ne+1),delta(2*ne+1),trial(2*ne+1)
    real(dp), allocatable :: jac(:,:)
    real(dp) :: grid(ne+1),left,right,xx,q,cn,exact,thermal,flux,flux_exact
    real(dp) :: alpha,base,trialnorm,scale,sgn,grad,etaq,oldcn,overshoot,backward
    integer :: e,a,b,ia,ib,row,col,it,bt,fld,idx(20),ns,ii
    logical :: accepted
    call clear_status(status)
    allocate(jac(2*ne+1,2*ne+1))
    call define_test_material(m)
    m%transport_supg=enabled;m%frozen_transport=.true.
    m%supg_directional_tau=.false.
    m%reference_n=1.0e20_dp;m%reference_p=1.0e20_dp
    m%mobility_n0=0.0_dp;m%mobility_p0=0.0_dp
    do a=1,3
      m%mobility_n0(a,a)=0.01_dp;m%mobility_p0(a,a)=0.01_dp
    end do
    thermal=BOLTZMANN_CONSTANT*m%reference_temperature/ELEMENTARY_CHARGE
    scale=ELEMENTARY_CHARGE*m%reference_n*0.01_dp*thermal
    call hex20_natural_coordinates(nat)
    fld=F_ETA_N;sgn=1.0_dp
    if(carrier==2) then
      fld=F_ETA_P;sgn=-1.0_dp
    end if
    do e=1,ne+1
      ! A priori boundary-layer mesh, mirrored for holes. No analytic initial guess.
      grid(e)=(real(e-1,dp)/ne)**3
      if(carrier==2) grid(e)=1.0_dp-(1.0_dp-real(e-1,dp)/ne)**3
    end do
    eta=log(1.5_dp);eta(1)=0.0_dp;eta(2*ne+1)=log(2.0_dp)
    load%joule_factor=0.0_dp
    do it=1,40
      trial=eta
      req%need_jacobian=.true.
      call assemble_bar()
      if(.not.status%is_ok()) return
      base=maxval(abs(res))
      if(base<1.0e-9_dp) exit
      call dense_solve(jac,-res,delta,status)
      if(.not.status%is_ok()) return
      alpha=min(1.0_dp,1.0_dp/max(maxval(abs(delta)),tiny(1.0_dp)))
      accepted=.false.
      do bt=1,24
        trial=eta+alpha*delta
        req%need_jacobian=.false.
        call assemble_bar()
        if(.not.status%is_ok()) return
        trialnorm=maxval(abs(res))
        if(trialnorm<base) then
          accepted=.true.;exit
        end if
        alpha=alpha*0.5_dp
      end do
      if(.not.accepted) exit
      eta=trial
    end do
    if(base>=1.0e-8_dp) then
      write(*,'(a,es12.4)') 'supg_bar_failed_residual=',base
      call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG bar Newton failed.');return
    end if
    error=0.0_dp;current_error=0.0_dp;oldcn=1.0_dp
    overshoot=0.0_dp;backward=0.0_dp
    do e=1,ne
      left=grid(e);right=grid(e+1)
      do ns=0,10
        q=-1.0_dp+0.2_dp*ns;xx=left+(q+1.0_dp)*0.5_dp*(right-left)
        etaq=0.5_dp*q*(q-1.0_dp)*eta(2*e-1)+(1.0_dp-q*q)*eta(2*e)+ &
          0.5_dp*q*(q+1.0_dp)*eta(2*e+1)
        grad=((q-0.5_dp)*eta(2*e-1)-2.0_dp*q*eta(2*e)+(q+0.5_dp)*eta(2*e+1))*2.0_dp/(right-left)
        cn=exp(etaq)
        if(carrier==1) then
          exact=1.0_dp+(1.0_dp-exp(-2.0_dp*pe*xx))/(1.0_dp-exp(-2.0_dp*pe))
          flux_exact=2.0_dp*pe*(2.0_dp-exp(-2.0_dp*pe))/(1.0_dp-exp(-2.0_dp*pe))
        else
          exact=1.0_dp+(exp(-2.0_dp*pe*(1.0_dp-xx))-exp(-2.0_dp*pe))/(1.0_dp-exp(-2.0_dp*pe))
          flux_exact=2.0_dp*pe*(1.0_dp-2.0_dp*exp(-2.0_dp*pe))/(1.0_dp-exp(-2.0_dp*pe))
        end if
        flux=cn*(2.0_dp*pe+sgn*grad)
        error=max(error,abs(cn-exact)/2.0_dp)
        current_error=max(current_error,abs(flux-flux_exact)/max(1.0_dp,abs(flux_exact)))
        overshoot=max(overshoot,1.0_dp-cn,cn-2.0_dp)
        backward=max(backward,oldcn-cn)
        oldcn=cn
      end do
    end do
    write(*,'(a,3(es12.4,a))') 'supg_bar_max_element_Pe=',pe*maxval(grid(2:)-grid(:ne)), &
      ',overshoot=',overshoot,',maximum_backward_step=',backward,','
    if(ne>=64 .and. (overshoot>1.0e-3_dp .or.backward>1.0e-3_dp)) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG resolved bar monotonicity tolerance failed.');return
    end if
    do a=1,20
      ii=1+nint(real(2*ne,dp)*(a-1)/19.0_dp)
      samples(a)=eta(ii)
    end do
  contains
    subroutine assemble_bar()
      res=0.0_dp;jac=0.0_dp
      do e=1,ne
        left=grid(e);right=grid(e+1)
        coords=0.5_dp*(nat+1.0_dp);coords(1,:)=left+coords(1,:)*(right-left)
        xe=0.0_dp
        do a=1,20
          idx(a)=2*e+nint(nat(1,a))
          xe(local_dof(a,F_PHI))=-2.0_dp*pe*thermal*coords(1,a)
          xe(local_dof(a,fld))=trial(idx(a))
        end do
        call evaluate_hex20_element(m,coords,xe,load,req,re,ke,diag,status)
        if(.not.status%is_ok()) return
        do a=1,20
          ia=idx(a);row=local_dof(a,fld)
          res(ia)=res(ia)+re(row)/scale
          if(req%need_jacobian) then
            do b=1,20
              ib=idx(b);col=local_dof(b,fld)
              jac(ia,ib)=jac(ia,ib)+ke(row,col)/scale
            end do
          end if
        end do
      end do
      res(1)=trial(1);res(2*ne+1)=trial(2*ne+1)-log(2.0_dp)
      jac(1,:)=0.0_dp;jac(1,1)=1.0_dp
      jac(2*ne+1,:)=0.0_dp;jac(2*ne+1,2*ne+1)=1.0_dp
    end subroutine assemble_bar
  end subroutine carrier_bar

  subroutine dense_solve(a,b,x,status)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(:)
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: c(:,:)
    real(dp) :: r(size(b)),tmp(size(b)),z,f
    integer :: n,i,j,p
    call clear_status(status)
    n=size(b);allocate(c(n,n));c=a;r=b
    do i=1,n
      p=i-1+maxloc(abs(c(i:n,i)),dim=1)
      if(abs(c(p,i))<tiny(1.0_dp)) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'SUPG test singular matrix.');return
      end if
      tmp=c(i,:);c(i,:)=c(p,:);c(p,:)=tmp
      z=r(i);r(i)=r(p);r(p)=z
      do j=i+1,n
        f=c(j,i)/c(i,i);c(j,i:n)=c(j,i:n)-f*c(i,i:n);r(j)=r(j)-f*r(i)
      end do
    end do
    x=0.0_dp
    do i=n,1,-1
      x(i)=(r(i)-dot_product(c(i,i+1:n),x(i+1:n)))/c(i,i)
    end do
  end subroutine dense_solve
end module supg_trace_iso_check_m

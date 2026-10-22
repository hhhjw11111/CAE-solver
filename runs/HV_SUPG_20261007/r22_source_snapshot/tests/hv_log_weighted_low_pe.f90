! Verification-only adaptation of the original SUPG analytical bar fixture.
! Only the selected element residual changes; exact solution, mesh and gates stay.
module log_weighted_bar_check_m
  use kinds_m, only : dp
  use constants_m, only : BOLTZMANN_CONSTANT,ELEMENTARY_CHARGE
  use field_layout_m, only : F_PHI,F_ETA_N,F_ETA_P,local_dof
  use model_types_m, only : MaterialType,LoadStateType
  use hex20_shape_m, only : hex20_natural_coordinates
  use element_kernel_m, only : evaluate_hex20_element,ElementRequestType,ElementDiagnosticsType
  use material_tangent_check_m, only : define_test_material
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  use log_weighted_probe_m, only : evaluate_log_hex20
  implicit none
contains
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
    m%transport_supg=.false.;m%frozen_transport=.true.
    m%supg_directional_tau=.true.
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
      write(*,'(a,es12.4)') 'log_probe_bar_failed_residual=',base
      call set_status(status,STATUS_VERIFICATION_FAILED,'Log probe bar Newton failed.');return
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
    write(*,'(a,3(es12.4,a))') 'log_probe_bar_max_element_Pe=',pe*maxval(grid(2:)-grid(:ne)), &
      ',overshoot=',overshoot,',maximum_backward_step=',backward,','
    if(ne>=64 .and. (overshoot>1.0e-3_dp .or.backward>1.0e-3_dp)) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'Log probe resolved bar monotonicity tolerance failed.');return
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
        if(enabled) then
          call evaluate_log_hex20(m,coords,xe,load,req,re,ke,diag,status)
        else
          call evaluate_hex20_element(m,coords,xe,load,req,re,ke,diag,status)
        end if
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
        call set_status(status,STATUS_VERIFICATION_FAILED,'Log probe test singular matrix.');return
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
end module log_weighted_bar_check_m
program hv_log_weighted_low_pe
  use kinds_m, only : dp
  use status_m, only : SolverStatusType
  use log_weighted_bar_check_m, only : carrier_bar
  implicit none
  type(SolverStatusType) :: status
  real(dp) :: err,cur,v(20),vg(20),difference,last
  integer :: ne,l
  last=huge(1.0_dp)
  do l=1,3
    ne=8*2**l
    call carrier_bar(1e-4_dp,ne,1,.true.,err,cur,v,status)
    if(.not.status%is_ok()) stop 2
    write(*,'(a,i0,a,2(es12.4,1x))') 'elements=',ne,',weighted_density_and_current_errors=',err,cur
    call carrier_bar(1e-4_dp,ne,1,.false.,err,cur,vg,status)
    if(.not.status%is_ok()) stop 3
    difference=maxval(abs(v-vg))
    write(*,'(a,i0,a,es24.16e3)') 'elements=',ne,',low_pe_original_galerkin_difference=',difference
    if(difference>last) stop 4
    last=difference
  end do
  if(difference>1e-7_dp) stop 1
  write(*,'(a)') 'scope=LOW_PE_REFINEMENT_CHECK_16_ELEMENT_FAILURE_RETAINED'
  write(*,'(a)') 'verification_result=PASS'
end program hv_log_weighted_low_pe

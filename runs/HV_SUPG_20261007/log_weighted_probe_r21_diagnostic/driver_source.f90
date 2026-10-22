program hv_log_weighted_probe
  use kinds_m, only : dp
  use constants_m, only : BOLTZMANN_CONSTANT,ELEMENTARY_CHARGE
  use model_types_m, only : MaterialType,LoadStateType
  use status_m, only : SolverStatusType
  use field_layout_m, only : F_PHI,F_ETA_N,F_ETA_P,F_THETA,local_dof
  use hex20_shape_m, only : hex20_natural_coordinates
  use material_tangent_check_m, only : define_test_material
  use element_kernel_m, only : ElementRequestType,ElementDiagnosticsType,evaluate_hex20_element
  use log_weighted_probe_m, only : evaluate_log_hex20
  implicit none
  type(MaterialType) :: m
  type(LoadStateType) :: load
  type(SolverStatusType) :: status
  type(ElementRequestType) :: request
  type(ElementDiagnosticsType) :: diagnostics
  real(dp) :: nat(3,20),x(3,20),u(140),d(140),r(140),k(140,140),rp(140),rm(140),unused(140,140)
  real(dp) :: original_r(140),original_k(140,140),action(140),fd(140),best(7),steps(5),scales(7)
  real(dp) :: err,den,worst,thermal,s,eta,normal_n,normal_p
  integer :: i,j,f,t,geometry,problem
  call define_test_material(m)
  m%reference_n=1.0e20_dp;m%reference_p=2.0e19_dp;m%reference_temperature=300.0_dp
  m%beta_n=2.0_dp;m%beta_p=2.3_dp
  m%transport_supg=.false.;m%linear_carriers=.false.;m%frozen_transport=.false.
  m%mobility_n0=0.0_dp;m%mobility_p0=0.0_dp
  do i=1,3
    m%mobility_n0(i,i)=0.01_dp*real(i,dp);m%mobility_p0(i,i)=0.008_dp*real(4-i,dp)
  end do
  load%voltage_factor=1.0_dp;load%joule_factor=1.0_dp
  call hex20_natural_coordinates(nat)
  scales=[1e-10_dp,1e-10_dp,1e-10_dp,0.1_dp,0.1_dp,0.1_dp,10.0_dp]
  steps=[1e-2_dp,1e-3_dp,1e-4_dp,1e-5_dp,1e-6_dp]
  worst=0.0_dp
  do geometry=1,2
    x=0.5e-6_dp*(nat+1.0_dp)
    if(geometry==2) then
      x(1,:)=x(1,:)+0.13e-6_dp*nat(2,:)*nat(3,:)
      x(2,:)=x(2,:)+0.08e-6_dp*nat(1,:)*nat(3,:)
    end if
    do i=1,140
      f=mod(i-1,7)+1
      u(i)=scales(f)*(0.3_dp*sin(real(i,dp)*0.23_dp)+0.1_dp*cos(real(i,dp)*0.4_dp))
    end do
    request%need_jacobian=.true.
    call evaluate_log_hex20(m,x,u,load,request,r,k,diagnostics,status)
    call require_ok()
    call evaluate_hex20_element(m,x,u,load,request,original_r,original_k,diagnostics,status)
    call require_ok()
    do f=1,7
      if(f==F_ETA_N.or.f==F_ETA_P) cycle
      if(any(r(f::7)/=original_r(f::7)).or.any(k(f::7,:)/=original_k(f::7,:))) stop 2
    end do
    ! Each input field separately exercises the full 7x7 row/column blocks.
    do j=1,7
      d=0.0_dp
      do i=1,20
        d(local_dof(i,j))=scales(j)*(sin(0.7_dp*i)+0.3_dp*cos(0.13_dp*i))
      end do
      action=matmul(k,d);best=huge(1.0_dp);request%need_jacobian=.false.
      do t=1,5
        call evaluate_log_hex20(m,x,u+steps(t)*d,load,request,rp,unused,diagnostics,status)
        call require_ok()
        call evaluate_log_hex20(m,x,u-steps(t)*d,load,request,rm,unused,diagnostics,status)
        call require_ok()
        fd=(rp-rm)/(2.0_dp*steps(t))
        do f=1,7
          den=max(maxval(abs(fd(f::7))),maxval(abs(action(f::7))),1e-200_dp)
          err=maxval(abs(fd(f::7)-action(f::7)))/den
          best(f)=min(best(f),err)
        end do
      end do
      worst=max(worst,maxval(best))
      write(*,'(a,i0,a,i0,a,es12.4)') 'geometry=',geometry,',input_field=',j,',maximum_block_error=',maxval(best)
      write(*,'(a,7(es12.4,1x))') 'row_block_errors=',best
      if(maxval(best)>5e-7_dp) stop 3
    end do
  end do
  thermal=BOLTZMANN_CONSTANT*300.0_dp/ELEMENTARY_CHARGE
  normal_n=ELEMENTARY_CHARGE*m%reference_n*maxval(m%mobility_n0)*thermal*1e-6_dp
  normal_p=ELEMENTARY_CHARGE*m%reference_p*maxval(m%mobility_p0)*thermal*1e-6_dp
  do problem=1,3
    x=0.5e-6_dp*(nat+1.0_dp);u=0.0_dp
    if(problem==1) then
      x(1,:)=x(1,:)+0.13e-6_dp*nat(2,:)*nat(3,:)
      x(2,:)=x(2,:)+0.08e-6_dp*nat(1,:)*nat(3,:)
    end if
    do i=1,20
      if(problem==1) then
        eta=1.1_dp*nat(1,i)+0.4_dp*nat(2,i)**2+0.3_dp*nat(2,i)*nat(3,i)
        u(local_dof(i,F_PHI))=thermal*eta
      else
        s=x(3,i)/1e-6_dp
        if(problem==3) s=(x(1,i)+0.35_dp*x(2,i)-0.2_dp*x(3,i))/1e-6_dp
        eta=3.0_dp*s
        u(local_dof(i,F_PHI))=BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE*3.0_dp*(300.0_dp*s+20.0_dp*s*s)
        u(local_dof(i,F_THETA))=40.0_dp*s
      end if
      u(local_dof(i,F_ETA_N))=eta;u(local_dof(i,F_ETA_P))=-eta
    end do
    call evaluate_log_hex20(m,x,u,load,request,r,unused,diagnostics,status)
    call require_ok()
    err=max(maxval(abs(r(F_ETA_N::7)))/normal_n,maxval(abs(r(F_ETA_P::7)))/normal_p)
    write(*,'(a,i0,a,es12.4)') 'exact_equilibrium=',problem,',normalized_carrier_residual=',err
    if(err>1e-10_dp) stop 4
  end do
  write(*,'(a,es24.16e3)') 'maximum_block_tangent_error=',worst
  write(*,'(a)') 'scope=VERIFICATION_ONLY_LOG_WEIGHTED_PROTOTYPE_NOT_PRODUCTION_OR_MESH_CERTIFICATE'
  write(*,'(a)') 'verification_result=PASS'
contains
  subroutine require_ok()
    if(.not.status%is_ok()) then
      write(*,'(a)') trim(status%message);stop 1
    end if
  end subroutine require_ok
end program hv_log_weighted_probe

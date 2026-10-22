! Exact carrier-current cancellations through the production HEX20 residual.
! These prescribed local fields need not solve the other five field equations.
program hv_transport_equilibrium
  use kinds_m, only : dp
  use constants_m, only : BOLTZMANN_CONSTANT,ELEMENTARY_CHARGE
  use field_layout_m, only : F_PHI,F_ETA_N,F_ETA_P,F_THETA,local_dof
  use model_types_m, only : MaterialType,LoadStateType
  use status_m, only : SolverStatusType
  use hex20_shape_m, only : hex20_natural_coordinates
  use material_tangent_check_m, only : define_test_material
  use element_kernel_m, only : ElementRequestType,ElementDiagnosticsType,evaluate_hex20_element
  implicit none
  type(MaterialType) :: material
  type(LoadStateType) :: load
  type(SolverStatusType) :: status
  type(ElementRequestType) :: request
  type(ElementDiagnosticsType) :: diagnostics
  real(dp) :: natural(3,20),coords(3,20),state(140),residual(140),jacobian(140,140)
  real(dp) :: thermal,s,eta,scale_n,scale_p,error_n,error_p,worst,coefficient(3)
  integer :: problem,node,enabled,directional,coefficient_id,i,checks
  real(dp), parameter :: length=1.0e-6_dp,limit=1.0e-10_dp
  call define_test_material(material)
  material%reference_temperature=300.0_dp;material%reference_n=1.0e20_dp;material%reference_p=2.0e19_dp
  material%beta_n=2.0_dp;material%beta_p=2.3_dp
  material%mobility_n0=0.0_dp;material%mobility_p0=0.0_dp
  do i=1,3
    material%mobility_n0(i,i)=0.01_dp*real(i,dp)
    material%mobility_p0(i,i)=0.008_dp*real(4-i,dp)
  end do
  material%linear_carriers=.false.;material%frozen_transport=.false.
  thermal=BOLTZMANN_CONSTANT*300.0_dp/ELEMENTARY_CHARGE
  coefficient=[0.5_dp,1.0_dp,2.0_dp]
  request%need_jacobian=.false.;load%joule_factor=0.0_dp
  call hex20_natural_coordinates(natural)
  worst=0.0_dp;checks=0
  do problem=1,3
    coords=length*0.5_dp*(natural+1.0_dp);state=0.0_dp
    if(problem==1) then
      ! Curved geometry, varying E/div(E), isothermal Boltzmann cancellation.
      coords(1,:)=coords(1,:)+length*0.13_dp*natural(2,:)*natural(3,:)
      coords(2,:)=coords(2,:)+length*0.08_dp*natural(1,:)*natural(3,:)
    end if
    do node=1,20
      if(problem==1) then
        eta=1.1_dp*natural(1,node)+0.4_dp*natural(2,node)**2+ &
            0.3_dp*natural(2,node)*natural(3,node)
        state(local_dof(node,F_PHI))=thermal*eta
      else
        ! s is affine, eta_n=3s, eta_p=-3s, T=300+40s, and
        ! phi=(kB/q)*3*(300s+20s^2), hence grad(phi)=(kBT/q)*grad(eta_n).
        s=coords(3,node)/length
        if(problem==3) s=(coords(1,node)+0.35_dp*coords(2,node)-0.2_dp*coords(3,node))/length
        eta=3.0_dp*s
        state(local_dof(node,F_THETA))=40.0_dp*s
        state(local_dof(node,F_PHI))=BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE* &
          3.0_dp*(300.0_dp*s+20.0_dp*s*s)
      end if
      state(local_dof(node,F_ETA_N))=eta;state(local_dof(node,F_ETA_P))=-eta
    end do
    scale_n=ELEMENTARY_CHARGE*material%reference_n*maxval(material%mobility_n0)*thermal*length* &
      exp(maxval(state(F_ETA_N::7)))
    scale_p=ELEMENTARY_CHARGE*material%reference_p*maxval(material%mobility_p0)*thermal*length* &
      exp(maxval(state(F_ETA_P::7)))
    do enabled=0,1
      material%transport_supg=enabled==1
      do directional=0,1
        material%supg_directional_tau=directional==1
        do coefficient_id=1,3
          material%supg_scale=coefficient(coefficient_id)
          call evaluate_hex20_element(material,coords,state,load,request,residual,jacobian,diagnostics,status)
          if(.not.status%is_ok()) then
            write(*,'(a)') trim(status%message);stop 1
          end if
          error_n=maxval(abs(residual(F_ETA_N::7)))/scale_n
          error_p=maxval(abs(residual(F_ETA_P::7)))/scale_p
          worst=max(worst,error_n,error_p);checks=checks+1
          write(*,'(a,i0,a,i0,a,i0,a,f3.1,a,2(es12.4,1x))') &
            'problem=',problem,',supg=',enabled,',directional=',directional, &
            ',coefficient=',material%supg_scale,',carrier_errors=',error_n,error_p
          if(max(error_n,error_p)>limit) stop 2
        end do
      end do
    end do
    ! Negative control: the same residual must detect a non-equilibrium field.
    state(local_dof(2,F_ETA_N))=state(local_dof(2,F_ETA_N))+0.01_dp
    state(local_dof(2,F_ETA_P))=state(local_dof(2,F_ETA_P))-0.01_dp
    call evaluate_hex20_element(material,coords,state,load,request,residual,jacobian,diagnostics,status)
    if(.not.status%is_ok()) stop 3
    if(min(maxval(abs(residual(F_ETA_N::7)))/scale_n,maxval(abs(residual(F_ETA_P::7)))/scale_p)<1.0e-7_dp) stop 4
  end do
  write(*,'(a,i0)') 'exact_cases=',checks
  write(*,'(a,es24.16e3)') 'maximum_normalized_carrier_residual=',worst
  write(*,'(a)') 'scope=LOCAL_CARRIER_RESIDUAL_CONSISTENCY_NOT_GLOBAL_MODEL_OR_MESH_CERTIFICATE'
  write(*,'(a)') 'verification_result=PASS'
end program hv_transport_equilibrium

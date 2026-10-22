program zno_point
 use kinds_m,only:dp
 use model_types_m,only:MaterialType,LoadStateType
 use status_m,only:SolverStatusType
 use zno_candidate_m,only:define_zno_candidate
 use material_point_m
 implicit none
 type(MaterialType)::m
 type(LoadStateType)::load
 type(MaterialPointStateType)::s,sp,sm
 type(MaterialPointResponseType)::rp,rm,rb
 type(MaterialPointTangentType)::t
 type(SolverStatusType)::status
 real(dp)::a(21,21),v(21),raw(21),ap(21),yp(21),ym(21),fd(21),err,h
 real(dp)::steps(9)
 integer::im,id,ih,i,u
 character(len=2)::models(3)
 models=['NT','N0','LQ']
 steps=[1e-2_dp,3e-3_dp,1e-3_dp,3e-4_dp,1e-4_dp,3e-5_dp,1e-5_dp,3e-6_dp,1e-6_dp]
 load%voltage_factor=1._dp;load%joule_factor=1._dp
 open(newunit=u,file='scan.csv',status='new')
 write(u,'(a)') 'model,direction,step,component,analytic,finite_difference,relative_error'
 do im=1,3
 call define_zno_candidate(models(im),m)
 call define_test_state(s)
 call evaluate_material_point(m,s,load,rb,t,status)
 if(.not.status%is_ok()) stop 1
 call material_tangent_matrix(t,a)
 write(*,'(a,a,a,es24.16,a,es24.16)') 'model=',models(im),',n=',rb%n,',p=',rb%p
 do id=1,3
 do i=1,21
 raw(i)=sin((0.37_dp+0.11_dp*id)*i)+cos((0.19_dp+0.07_dp*id)*i)
 end do
 raw=raw/sqrt(sum(raw*raw));v=raw*input_scales();ap=matmul(a,v)
 do ih=1,9
 h=steps(ih)
 call perturb_state(s,v,h,sp);call perturb_state(s,v,-h,sm)
 call evaluate_material_point(m,sp,load,rp,t,status)
 if(.not.status%is_ok()) stop 2
 call evaluate_material_point(m,sm,load,rm,t,status)
 if(.not.status%is_ok()) stop 3
 call pack_material_response(rp,yp);call pack_material_response(rm,ym)
 fd=(yp-ym)/(2*h)
 do i=1,21
 err=abs(fd(i)-ap(i))/max(abs(fd(i)),abs(ap(i)),1e-280_dp)
 write(u,'(a,a,i0,a,es24.16,a,i0,3(a,es24.16))') models(im),',',id,',',h,',',i,',',ap(i),',',fd(i),',',err
 end do
 end do
 end do
 end do
 close(u)
contains
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

end program zno_point

program compare_kernel
  use kinds_m, only: dp
  use model_types_m, only: MaterialType,LoadStateType
  use status_m, only: SolverStatusType
  use field_layout_m, only: local_dof
  use hex20_shape_m, only: hex20_natural_coordinates
  use material_tangent_check_m, only: define_test_material
  use element_kernel_m, only: ElementRequestType,ElementDiagnosticsType,evaluate_hex20_element
  implicit none
  type(MaterialType) :: m
  type(LoadStateType) :: load
  type(SolverStatusType) :: status
  type(ElementRequestType) :: request
  type(ElementDiagnosticsType) :: diagnostics
  real(dp) :: coords(3,20),nat(3,20),state(140),r(140),k(140,140),x,y,z,t0,t1,total
  integer :: case_id,node,rep,unit_id
  character(len=32) :: filename
  total=0
  do case_id=1,4
    call define_test_material(m)
    m%reference_p=2e20_dp;m%transport_supg=.true.
    m%supg_directional_tau=case_id/=1;m%frozen_transport=case_id==3
    call hex20_natural_coordinates(nat)
    coords=.5_dp*(nat+1)
    do node=1,20
      x=coords(1,node);y=coords(2,node);z=coords(3,node)
      coords(:,node)=coords(:,node)+[.02_dp*y*z,.01_dp*x*z,.015_dp*x*y]
    end do
    state=0
    do node=1,20
      x=coords(1,node);y=coords(2,node);z=coords(3,node)
      state(local_dof(node,1))=1e-4_dp*x+.2e-4_dp*y-.1e-4_dp*z
      state(local_dof(node,2))=-.3e-4_dp*x+.8e-4_dp*y+.2e-4_dp*z
      state(local_dof(node,3))=.1e-4_dp*x-.2e-4_dp*y+.6e-4_dp*z
      state(local_dof(node,4))=-1.2e5_dp*x+.8e5_dp*y-.6e5_dp*z+2e4_dp*x*y
      if(case_id==4) state(local_dof(node,4))=0
      state(local_dof(node,5))=.12_dp+.05_dp*x-.03_dp*y+.02_dp*z+.01_dp*x*z
      state(local_dof(node,6))=-.18_dp-.04_dp*x+.06_dp*y-.01_dp*z-.02_dp*y*z
      state(local_dof(node,7))=45+3*x-2*y+1.5_dp*z+2*x*y
    end do
    load%voltage_factor=1;load%joule_factor=.81_dp;request%need_jacobian=.true.
    call cpu_time(t0)
    do rep=1,40
      call evaluate_hex20_element(m,coords,state,load,request,r,k,diagnostics,status)
      if(.not.status%is_ok()) then
        print *,trim(status%message)
        stop 1
      end if
    end do
    call cpu_time(t1)
    total=total+t1-t0
    write(filename,'(a,i0,a)') 'case_',case_id,'.bin'
    open(newunit=unit_id,file=trim(filename),access='stream',form='unformatted',status='new')
    write(unit_id) r,k
    close(unit_id)
    write(*,'(a,i0,a,es16.8)') 'case=',case_id,',cpu_seconds=',t1-t0
  end do
  write(*,'(a,es16.8)') 'total_cpu_seconds=',total
end program compare_kernel

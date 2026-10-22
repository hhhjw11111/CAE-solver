module a2_hole_check_m
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE,F_PHI,F_ETA_N,F_ETA_P,F_THETA,global_dof
  use model_types_m, only : MaterialType,StateType,LoadStateType,NonlinearOptionsType
  use mesh_types_m, only : MeshType
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  use hex20_shape_m, only : hex20_natural_coordinates,hex20_shape,hex20_geometry
  use mesh_checks_m, only : validate_mesh
  use sparse_pattern_m, only : SparseAssemblyContextType,initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use newton_solver_m, only : NewtonReportType,solve_single_load_newton
  use independent_flux_audit_m, only : independent_gauss_rule
  implicit none
  private
  real(dp), parameter :: HALF_WIDTH=10.0e-6_dp,RADIUS=3.0e-6_dp,THICKNESS=1.0e-6_dp,E0=1.0e4_dp
  public :: verify_a2_hole
contains
  subroutine verify_a2_hole(directory,status)
    character(len=*), intent(in) :: directory
    type(SolverStatusType), intent(out) :: status
    integer, parameter :: nt(3)=[16,24,32],nr(3)=[2,3,4]
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(NonlinearOptionsType) :: options
    type(NewtonReportType) :: report
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: bc(:)
    real(dp) :: errors(5,3),minimum_det,volume,scales(7)
    integer :: level,i,nnode
    logical :: passed
    character(len=16) :: label
    call clear_status(status)
    write(*,'(a)') 'verification_id=A2_F1_insulated_hole'
    write(*,'(a)') 'model=production_nonlinear_two_carrier_isothermal_synthetic'
    materials(1)=MaterialType()
    materials(1)%name='A2_SYNTHETIC_ISOTHERMAL'
    materials(1)%reference_n=1.0e21_dp;materials(1)%reference_p=1.0e21_dp
    materials(1)%intrinsic_concentration=1.0e21_dp
    do i=1,6
      materials(1)%elastic(i,i)=1.0e11_dp
    end do
    do i=1,3
      materials(1)%permittivity(i,i)=1.0e-10_dp
      materials(1)%thermal_conductivity(i,i)=100.0_dp
      materials(1)%mobility_n0(i,i)=0.02_dp
      materials(1)%mobility_p0(i,i)=0.005_dp
    end do
    scales=[1.0e-10_dp,1.0e-10_dp,1.0e-10_dp,E0*HALF_WIDTH,0.1_dp,0.1_dp,1.0_dp]
    options%max_newton=30;options%residual_tolerance=1.0e-11_dp;options%update_tolerance=1.0e-11_dp
    load%voltage_factor=1.0_dp;load%joule_factor=1.0_dp
    passed=.true.
    do level=1,3
      write(label,'("L",i0)') level-1
      call build_hole(nt(level),nr(level),mesh)
      call validate_mesh(mesh,minimum_det,volume,status)
      if (.not. status%is_ok()) return
      call define_bc(mesh,bc)
      call initialize_sparse_context(mesh,context,status)
      if (.not. status%is_ok()) return
      nnode=mesh%number_of_nodes()
      if (allocated(state%x)) deallocate(state%x)
      allocate(state%x(NDOF_NODE*nnode));state%x=0.0_dp
      call solve_single_load_newton(mesh,materials,load,context,bc,scales,options,state,report,status)
      if (.not. status%is_ok()) return
      call measure(mesh,state,nr(level),trim(directory),trim(label),errors(:,level),status)
      if (.not. status%is_ok()) return
      call metric(label,'phi_relative_L2',errors(1,level))
      call metric(label,'E_relative_L2',errors(2,level))
      call metric(label,'wall_normal_E_RMS_over_E0',errors(3,level))
      call metric(label,'wall_Et_peak_over_E0',errors(4,level))
      call metric(label,'wall_Q_peak_over_Qinf',errors(5,level))
      call metric(label,'max_relative_density_shift',maxval(abs(exp(state%x(F_ETA_N::NDOF_NODE))-1.0_dp)))
      call metric(label,'residual',report%final_residual_rms)
      call metric(label,'update',report%final_update_rms)
      write(*,'(a,i0)') 'a2_'//trim(label)//'_nodes=',nnode
      write(*,'(a,i0)') 'a2_'//trim(label)//'_newton_iterations=',report%iterations
      passed=passed .and. report%converged .and. report%maximum_perturbed_pivots==0
    end do
    passed=passed .and. errors(1,3)<1.0e-3_dp .and. errors(2,3)<0.02_dp .and. &
      errors(3,3)<0.01_dp .and. abs(errors(4,3)/2.0_dp-1.0_dp)<0.02_dp .and. &
      abs(errors(5,3)/4.0_dp-1.0_dp)<0.02_dp
    passed=passed .and. all(errors(1:3,2)<errors(1:3,1)) .and. all(errors(1:3,3)<errors(1:3,2))
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'A2 frozen analytic criteria not satisfied.')
    end if
  end subroutine verify_a2_hole

  subroutine build_hole(nt,nr,mesh)
    integer, intent(in) :: nt,nr
    type(MeshType), intent(out) :: mesh
    integer, allocatable :: ids(:,:,:)
    real(dp), allocatable :: coordinates(:,:)
    real(dp) :: natural(3,20),angle,outer,rho,r,pi
    integer :: er,et,j,ir,it,iz,node,element
    pi=acos(-1.0_dp)
    call hex20_natural_coordinates(natural)
    allocate(ids(0:2*nr,0:2*nt-1,0:2),coordinates(3,20*nr*nt))
    allocate(mesh%connectivity(20,nr*nt),mesh%material_id(nr*nt))
    ids=0;node=0;element=0;mesh%material_id=1
    do et=0,nt-1
      do er=0,nr-1
        element=element+1
        do j=1,20
          ir=2*er+nint(natural(1,j))+1
          it=modulo(2*et+nint(natural(2,j))+1,2*nt)
          iz=nint(natural(3,j))+1
          if (ids(ir,it,iz)==0) then
            node=node+1;ids(ir,it,iz)=node
            angle=pi*real(it,dp)/real(nt,dp)
            outer=HALF_WIDTH/max(abs(cos(angle)),abs(sin(angle)))
            rho=(real(ir,dp)/real(2*nr,dp))**2
            ! Radial midside nodes are physical midpoints, not re-graded points.
            if (modulo(ir,2)==1) rho=0.5_dp*((real(ir-1,dp)/real(2*nr,dp))**2+ &
              (real(ir+1,dp)/real(2*nr,dp))**2)
            r=RADIUS+rho*(outer-RADIUS)
            coordinates(:,node)=[r*cos(angle),r*sin(angle),THICKNESS*0.5_dp*real(iz-1,dp)]
          end if
          mesh%connectivity(j,element)=ids(ir,it,iz)
        end do
      end do
    end do
    mesh%coordinates=coordinates(:,1:node)
  end subroutine build_hole

  subroutine define_bc(mesh,bc)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), allocatable, intent(out) :: bc(:)
    integer, allocatable :: nodes(:),outer(:)
    real(dp) :: phi,efield(3)
    integer :: i,j
    nodes=[(i,i=1,mesh%number_of_nodes())]
    outer=pack(nodes,abs(maxval(abs(mesh%coordinates(1:2,:)),dim=1)-HALF_WIDTH)<1.0e-15_dp)
    allocate(bc(6+size(outer)))
    do i=1,3
      bc(i)%field_id=i;bc(i)%node_ids=nodes
    end do
    bc(4)%field_id=F_THETA;bc(4)%node_ids=nodes
    bc(5)%field_id=F_ETA_N;bc(5)%node_ids=outer
    bc(6)%field_id=F_ETA_P;bc(6)%node_ids=outer
    do j=1,size(outer)
      call exact(mesh%coordinates(:,outer(j)),phi,efield)
      bc(6+j)%field_id=F_PHI;bc(6+j)%node_ids=[outer(j)];bc(6+j)%value_end=phi
    end do
  end subroutine define_bc

  subroutine exact(x,phi,efield)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: phi,efield(3)
    real(dp) :: r2
    r2=sum(x(1:2)**2)
    phi=-E0*x(1)*(1.0_dp+RADIUS**2/r2)
    efield=E0*[1.0_dp+RADIUS**2*(x(2)**2-x(1)**2)/r2**2, &
      -2.0_dp*RADIUS**2*x(1)*x(2)/r2**2,0.0_dp]
  end subroutine exact

  subroutine measure(mesh,state,nr,directory,label,errors,status)
    type(MeshType), intent(in) :: mesh
    type(StateType), intent(in) :: state
    integer, intent(in) :: nr
    character(len=*), intent(in) :: directory,label
    real(dp), intent(out) :: errors(5)
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: gp(:),gw(:)
    real(dp) :: coords(3,20),shape(20),gradient(3,20),dn(3,20),jac(3,3),phi_nodes(20)
    real(dp) :: x(3),phi_ex,e_ex(3),e_num(3),det,w,area,normal(3),en,et,qratio,eta,t
    real(dp) :: sums(6),peak_et,peak_q,phi_num
    integer :: element,i,j,k,node,unit_nodes,unit_wall,io
    logical :: saving
    call independent_gauss_rule(5,gp,gw,status)
    if (.not. status%is_ok()) return
    saving=len_trim(directory)>0
    if (saving) then
      open(newunit=unit_nodes,file=trim(directory)//'/A2_'//trim(label)//'_nodes.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create A2 nodal file.');return
      end if
      write(unit_nodes,'(a)') 'x_m,y_m,z_m,phi_V,phi_exact_V,eta_n,eta_p'
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(:,node);call exact(x,phi_ex,e_ex)
        write(unit_nodes,'(*(es24.16e3,:,","))') x,state%x(global_dof(node,F_PHI)),phi_ex, &
          state%x(global_dof(node,F_ETA_N)),state%x(global_dof(node,F_ETA_P))
      end do
      close(unit_nodes)
      open(newunit=unit_wall,file=trim(directory)//'/A2_'//trim(label)//'_wall.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create A2 wall file.');return
      end if
      write(unit_wall,'(a)') 'angle_rad,x_m,y_m,En_over_E0,Et_over_E0,Q_over_Qinf,Et_exact_over_E0,Q_exact_over_Qinf'
    end if
    sums=0.0_dp;peak_et=0.0_dp;peak_q=0.0_dp
    do element=1,mesh%number_of_elements()
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      do node=1,20
        phi_nodes(node)=state%x(global_dof(mesh%connectivity(node,element),F_PHI))
      end do
      do k=1,5
        do j=1,5
          do i=1,5
            call hex20_geometry(coords,gp(i),gp(j),gp(k),shape,gradient,det,status)
            if (.not. status%is_ok()) return
            w=gw(i)*gw(j)*gw(k)*det;x=matmul(coords,shape)
            call exact(x,phi_ex,e_ex)
            phi_num=dot_product(shape,phi_nodes);e_num=-matmul(gradient,phi_nodes)
            sums(1:4)=sums(1:4)+w*[(phi_num-phi_ex)**2,phi_ex**2,sum((e_num-e_ex)**2),sum(e_ex**2)]
          end do
        end do
      end do
      if (modulo(element-1,nr)/=0) cycle
      do k=1,5
        do j=1,5
          call wall_point(gp(j),gp(k))
          if (.not. status%is_ok()) return
          w=gw(j)*gw(k)*area
          sums(5)=sums(5)+w*(en/E0)**2;sums(6)=sums(6)+w
        end do
      end do
      ! Include end points so the fixed 90/270-degree maxima are sampled.
      do j=0,20
        t=-1.0_dp+0.1_dp*real(j,dp)
        call wall_point(t,0.0_dp)
        if (.not. status%is_ok()) return
        peak_et=max(peak_et,et/E0);peak_q=max(peak_q,qratio)
        eta=atan2(x(2),x(1))
        if (saving) write(unit_wall,'(*(es24.16e3,:,","))') eta,x(1:2),en/E0,et/E0,qratio, &
          2.0_dp*abs(sin(eta)),4.0_dp*sin(eta)**2
      end do
    end do
    if (saving) close(unit_wall)
    errors=[sqrt(sums(1)/sums(2)),sqrt(sums(3)/sums(4)),sqrt(sums(5)/sums(6)),peak_et,peak_q]
  contains
    subroutine wall_point(t,z)
      real(dp), intent(in) :: t,z
      call hex20_geometry(coords,-1.0_dp,t,z,shape,gradient,det,status)
      if (.not. status%is_ok()) return
      call hex20_shape(-1.0_dp,t,z,shape,dn)
      jac=matmul(coords,transpose(dn))
      normal=[jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3), &
        jac(3,2)*jac(1,3)-jac(1,2)*jac(3,3),jac(1,2)*jac(2,3)-jac(2,2)*jac(1,3)]
      area=sqrt(sum(normal**2));normal=-normal/area
      x=matmul(coords,shape);e_num=-matmul(gradient,phi_nodes)
      en=dot_product(e_num,normal);et=sqrt(max(0.0_dp,sum(e_num**2)-en**2))
      qratio=sum(e_num**2)/E0**2
    end subroutine wall_point
  end subroutine measure

  subroutine metric(label,name,value)
    character(len=*), intent(in) :: label,name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') 'a2_'//trim(label)//'_'//trim(name)//'=',value
  end subroutine metric
end module a2_hole_check_m

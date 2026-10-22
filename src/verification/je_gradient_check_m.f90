module je_gradient_check_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE, BOLTZMANN_CONSTANT
  use field_layout_m, only : NDOF_NODE, F_PHI, F_ETA_N, F_ETA_P, F_THETA, global_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use model_types_m, only : MaterialType, StateType, LoadStateType, NonlinearOptionsType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use hex20_shape_m, only : hex20_natural_coordinates, hex20_geometry
  use mesh_checks_m, only : validate_mesh
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use newton_solver_m, only : NewtonReportType, solve_single_load_newton
  use independent_flux_audit_m, only : IndependentFluxReportType, &
    integrate_independent_fluxes, independent_gauss_rule
  use material_point_m, only : MaterialPointStateType, MaterialPointResponseType, &
    MaterialPointTangentType, evaluate_material_point
  implicit none
  private
  real(dp), parameter :: LENGTH=20.0e-6_dp, WIDTH=20.0e-6_dp, THICKNESS=1.0e-6_dp
  real(dp), parameter :: DENSITY=1.0e21_dp, MU_N=0.05_dp, MU_P=0.02_dp, KAPPA=10.0_dp
  real(dp), parameter :: T0=300.0_dp, VOLTAGE=0.02_dp, ETA_SCALE=0.05_dp
  real(dp), parameter :: VT=BOLTZMANN_CONSTANT*T0/ELEMENTARY_CHARGE
  real(dp), parameter :: FIELD_TOL=1.0e-8_dp, POWER_TOL=1.0e-8_dp, QUAD_TOL=1.0e-10_dp
  real(dp), parameter :: FIELD_SCALES(7)=[1.0e-10_dp,1.0e-10_dp,1.0e-10_dp, &
    VOLTAGE,ETA_SCALE,ETA_SCALE,1.0e-4_dp]
  integer, parameter :: COUNTS(3)=[1,2,4]
  public :: verify_je_gradient

contains

  subroutine verify_je_gradient(status)
    type(SolverStatusType), intent(out) :: status
    integer :: sign_index,level,total
    real(dp) :: slope
    logical :: passed,one_pass
    character(len=48) :: label
    call clear_status(status)
    passed=.true.; total=0
    write(*,'(a)') 'verification_id=JEGRAD_production_nonzero_diffusion_work'
    write(*,'(a)') 'jegrad_scope=synthetic_LQ_frozen_transport_exact_solution_not_ZnO_validation'
    write(*,'(a)') 'jegrad_interpolation=linear_phi_eta_quadratic_theta_exactly_representable'
    write(*,'(a)') 'jegrad_LQ_density=reference_density_times_one_plus_eta'
    write(*,'(a)') 'jegrad_diffusion_work=independent_reference_density_decomposition'
    call metric('protocol','field_tolerance',FIELD_TOL)
    call metric('protocol','power_tolerance',POWER_TOL)
    call metric('protocol','quadrature_tolerance',QUAD_TOL)
    do sign_index=-1,1
      slope=real(sign_index,dp)*ETA_SCALE
      do level=1,size(COUNTS)
        write(label,'(a,i0,a,i0)') 's',sign_index,'_nx',COUNTS(level)
        call solve_case(trim(label),COUNTS(level),1,1,slope,one_pass,status)
        if (.not. status%is_ok()) return
        passed=passed .and. one_pass; total=total+1
      end do
    end do
    call solve_case('transverse_2x2x2',2,2,2,ETA_SCALE,one_pass,status)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass; total=total+1
    call verify_nt_material(one_pass,status)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass
    write(*,'(a,i0)') 'jegrad_solved_cases=',total
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'JEGRAD exceeded a frozen field/power/tangent criterion.')
    end if
  end subroutine verify_je_gradient

  subroutine solve_case(label,nx,ny,nz,slope,passed,status)
    character(len=*), intent(in) :: label
    integer, intent(in) :: nx,ny,nz
    real(dp), intent(in) :: slope
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1),old_materials(1)
    type(StateType) :: state,bad_state
    type(LoadStateType) :: load
    type(NonlinearOptionsType) :: options
    type(NewtonReportType) :: newton
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: boundaries(:)
    type(IndependentFluxReportType) :: flux,flux4,old_flux,bad_flux
    real(dp) :: jn,jp,qje,qdrift,qdiff,power,power_scale,area,volume,min_det
    real(dp) :: errors(8),point_errors(3),flux_errors(11),quad_error,geometry_error
    real(dp) :: numerical_drift,negative_heat,negative_source,diffusion_error,peak
    integer :: expected_exterior,expected_internal
    call clear_status(status)
    passed=.false.
    call build_box(nx,ny,nz,mesh)
    call define_material(materials(1))
    call analytic_constants(slope,jn,jp,qje,qdrift,qdiff)
    area=WIDTH*THICKNESS; power=qje*LENGTH*area; power_scale=qdrift*LENGTH*area
    peak=qje*LENGTH**2/(8.0_dp*KAPPA)
    call validate_mesh(mesh,min_det,volume,status)
    if (.not. status%is_ok()) return
    call define_boundaries(mesh,slope,boundaries)
    call initialize_sparse_context(mesh,context,status)
    if (.not. status%is_ok()) return
    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()))
    state%x=0.0_dp
    load%voltage_factor=1.0_dp; load%joule_factor=1.0_dp
    options%max_newton=30
    options%residual_tolerance=1.0e-12_dp; options%update_tolerance=1.0e-12_dp
    call solve_single_load_newton(mesh,materials,load,context,boundaries,FIELD_SCALES, &
      options,state,newton,status)
    if (.not. status%is_ok()) then
      write(*,'(a)') 'jegrad_failed_case='//trim(label)
      call metric(label,'failed_residual',newton%final_residual_rms)
      return
    end if
    call field_errors(mesh,materials(1),state,slope,errors,point_errors,numerical_drift,status)
    if (.not. status%is_ok()) return
    call integrate_independent_fluxes(mesh,materials,state,load,5,flux,status)
    if (.not. status%is_ok()) return
    call integrate_independent_fluxes(mesh,materials,state,load,4,flux4,status)
    if (.not. status%is_ok()) return
    flux_errors(1)=maxval(abs(flux%electron_out-[-jn*area,jn*area,0.0_dp]))/(abs(jn)*area)
    flux_errors(2)=maxval(abs(flux%hole_out-[-jp*area,jp*area,0.0_dp]))/(abs(jp)*area)
    flux_errors(3)=maxval(abs(flux%heat_out-[power/2.0_dp,power/2.0_dp,0.0_dp]))/power_scale
    flux_errors(4)=abs(flux%joule-power)/power_scale
    flux_errors(5)=abs(flux%terminal_power-power)/power_scale
    flux_errors(6)=abs(flux%total_current_work-power)/power_scale
    flux_errors(7)=abs(sum(flux%heat_out)-flux%joule)/power_scale
    flux_errors(8)=abs(flux%terminal_power-flux%total_current_work)/power_scale
    flux_errors(9)=abs(sum(flux%electron_out))/(abs(jn)*area)
    flux_errors(10)=abs(sum(flux%hole_out))/(abs(jp)*area)
    flux_errors(11)=abs(flux%charge)/(ELEMENTARY_CHARGE*DENSITY*LENGTH*area)
    diffusion_error=abs(flux%joule-numerical_drift-qdiff*LENGTH*area)/power_scale
    quad_error=max(abs(flux%joule-flux4%joule)/power_scale, &
      abs(flux%terminal_power-flux4%terminal_power)/power_scale, &
      maxval(abs(flux%heat_out-flux4%heat_out))/power_scale, &
      maxval(abs(flux%electron_out-flux4%electron_out))/(abs(jn)*area), &
      maxval(abs(flux%hole_out-flux4%hole_out))/(abs(jp)*area))
    geometry_error=max(abs(flux%volume/(LENGTH*area)-1.0_dp), &
      abs(flux%area/(2.0_dp*(LENGTH*WIDTH+LENGTH*THICKNESS+area))-1.0_dp), &
      sqrt(sum(flux%area_normal**2))/flux%area)
    expected_exterior=2*(nx*ny+ny*nz+nz*nx)
    expected_internal=(nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1)
    passed=newton%converged .and. newton%maximum_perturbed_pivots == 0 .and. &
      maxval(errors) < FIELD_TOL .and. maxval(point_errors) < FIELD_TOL .and. &
      maxval(flux_errors) < POWER_TOL .and. diffusion_error < POWER_TOL .and. &
      quad_error < QUAD_TOL .and. geometry_error < QUAD_TOL .and. &
      flux%exterior_faces == expected_exterior .and. flux%internal_faces == expected_internal
    call metric(label,'slope',slope)
    call metric(label,'phi_L2_scaled',errors(1))
    call metric(label,'E_L2_scaled',errors(2))
    call metric(label,'eta_n_L2_scaled',errors(3))
    call metric(label,'eta_p_L2_scaled',errors(4))
    call metric(label,'grad_eta_n_L2_scaled',errors(5))
    call metric(label,'grad_eta_p_L2_scaled',errors(6))
    call metric(label,'theta_L2_relative',errors(7))
    call metric(label,'grad_theta_L2_relative',errors(8))
    call metric(label,'point_Jn_relative',point_errors(1))
    call metric(label,'point_Jp_relative',point_errors(2))
    call metric(label,'point_Q_relative',point_errors(3))
    call metric(label,'max_flux_error',maxval(flux_errors))
    call metric(label,'diffusion_work_error',diffusion_error)
    call metric(label,'quadrature_4_vs_5',quad_error)
    call metric(label,'geometry_error',geometry_error)
    call metric(label,'joule_W',flux%joule)
    call metric(label,'terminal_power_W',flux%terminal_power)
    call metric(label,'diffusion_work_W',flux%joule-numerical_drift)
    call metric(label,'analytic_diffusion_work_W',qdiff*LENGTH*area)
    call metric(label,'analytic_peak_K',peak)
    call metric(label,'numerical_peak_K',maxval(state%x(F_THETA::NDOF_NODE)))
    call metric(label,'newton_residual',newton%final_residual_rms)
    write(*,'(a,i0)') 'jegrad_'//trim(label)//'_newton_iterations=',newton%iterations
    if (nx == 2 .and. ny == 1 .and. slope > 0.0_dp) then
      old_materials=materials; old_materials(1)%heat_source_current_work=.false.
      call integrate_independent_fluxes(mesh,old_materials,state,load,5,old_flux,status)
      if (.not. status%is_ok()) return
      negative_source=abs(old_flux%joule-power)/power_scale
      allocate(bad_state%x(size(state%x))); bad_state%x=state%x
      bad_state%x(F_THETA::NDOF_NODE)=1.1_dp*bad_state%x(F_THETA::NDOF_NODE)
      call integrate_independent_fluxes(mesh,materials,bad_state,load,5,bad_flux,status)
      if (.not. status%is_ok()) return
      negative_heat=abs(sum(bad_flux%heat_out)-bad_flux%joule)/power_scale
      passed=passed .and. negative_source > 0.01_dp .and. negative_heat > 0.05_dp
      call metric(label,'negative_old_source_defect',negative_source)
      call metric(label,'negative_heat_control_defect',negative_heat)
    end if
    write(*,'(a,l1)') 'jegrad_'//trim(label)//'_passed=',passed
  end subroutine solve_case

  subroutine analytic_constants(slope,jn,jp,qje,qdrift,qdiff)
    real(dp), intent(in) :: slope
    real(dp), intent(out) :: jn,jp,qje,qdrift,qdiff
    real(dp) :: efield
    efield=VOLTAGE/LENGTH
    jn=ELEMENTARY_CHARGE*DENSITY*MU_N*(efield+VT*slope/LENGTH)
    jp=ELEMENTARY_CHARGE*DENSITY*MU_P*(efield-VT*slope/LENGTH)
    qdrift=ELEMENTARY_CHARGE*DENSITY*(MU_N+MU_P)*efield**2
    qdiff=ELEMENTARY_CHARGE*DENSITY*(MU_N-MU_P)*VT*efield*slope/LENGTH
    qje=qdrift+qdiff
  end subroutine analytic_constants

  subroutine define_material(material)
    type(MaterialType), intent(out) :: material
    integer :: i
    material=MaterialType()
    material%name='JEGRAD_SYNTHETIC_LQ_NEUTRAL'
    material%reference_temperature=T0
    material%reference_n=DENSITY; material%reference_p=DENSITY
    material%intrinsic_concentration=DENSITY
    material%frozen_transport=.true.; material%linear_carriers=.true.
    material%heat_source_current_work=.true.
    do i=1,6
      material%elastic(i,i)=1.0e11_dp
    end do
    do i=1,3
      material%permittivity(i,i)=8.0e-11_dp
      material%thermal_conductivity(i,i)=KAPPA
      material%mobility_n0(i,i)=MU_N; material%mobility_p0(i,i)=MU_P
    end do
  end subroutine define_material

  subroutine build_box(nx,ny,nz,mesh)
    integer, intent(in) :: nx,ny,nz
    type(MeshType), intent(out) :: mesh
    integer, allocatable :: ids(:,:,:)
    real(dp), allocatable :: coords(:,:)
    real(dp) :: natural(3,HEX20_NNODE)
    integer :: ix,iy,iz,node,gx,gy,gz,number,element
    call hex20_natural_coordinates(natural)
    allocate(ids(0:2*nx,0:2*ny,0:2*nz),coords(3,20*nx*ny*nz))
    allocate(mesh%connectivity(20,nx*ny*nz),mesh%material_id(nx*ny*nz))
    ids=0; number=0; element=0; mesh%material_id=1
    do iz=0,nz-1
      do iy=0,ny-1
        do ix=0,nx-1
          element=element+1
          do node=1,20
            gx=2*ix+nint(natural(1,node))+1
            gy=2*iy+nint(natural(2,node))+1
            gz=2*iz+nint(natural(3,node))+1
            if (ids(gx,gy,gz) == 0) then
              number=number+1; ids(gx,gy,gz)=number
              coords(:,number)=[LENGTH*real(gx,dp)/real(2*nx,dp), &
                WIDTH*(real(gy,dp)/real(2*ny,dp)-0.5_dp), &
                THICKNESS*(real(gz,dp)/real(2*nz,dp)-0.5_dp)]
            end if
            mesh%connectivity(node,element)=ids(gx,gy,gz)
          end do
        end do
      end do
    end do
    mesh%coordinates=coords(:,1:number)
  end subroutine build_box

  subroutine define_boundaries(mesh,slope,boundaries)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: slope
    type(DirichletBCType), allocatable, intent(out) :: boundaries(:)
    integer, allocatable :: nodes(:),left(:),right(:)
    integer :: i
    nodes=[(i,i=1,mesh%number_of_nodes())]
    left=pack(nodes,abs(mesh%coordinates(1,:)) < 1.0e-16_dp)
    right=pack(nodes,abs(mesh%coordinates(1,:)-LENGTH) < 1.0e-16_dp)
    allocate(boundaries(11))
    do i=1,3
      call set_boundary(boundaries(i),i,nodes,0.0_dp)
    end do
    call set_boundary(boundaries(4),F_PHI,left,0.0_dp)
    call set_boundary(boundaries(5),F_PHI,right,-VOLTAGE)
    call set_boundary(boundaries(6),F_ETA_N,left,0.0_dp)
    call set_boundary(boundaries(7),F_ETA_N,right,slope)
    call set_boundary(boundaries(8),F_ETA_P,left,0.0_dp)
    call set_boundary(boundaries(9),F_ETA_P,right,slope)
    call set_boundary(boundaries(10),F_THETA,left,0.0_dp)
    call set_boundary(boundaries(11),F_THETA,right,0.0_dp)
  end subroutine define_boundaries

  subroutine set_boundary(boundary,field,nodes,value)
    type(DirichletBCType), intent(out) :: boundary
    integer, intent(in) :: field,nodes(:)
    real(dp), intent(in) :: value
    boundary%field_id=field; boundary%node_ids=nodes
    boundary%value_start=0.0_dp; boundary%value_end=value
  end subroutine set_boundary

  subroutine field_errors(mesh,material,state,slope,errors,point_errors,drift_integral,status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    real(dp), intent(in) :: slope
    real(dp), intent(out) :: errors(8),point_errors(3),drift_integral
    type(SolverStatusType), intent(out) :: status
    type(MaterialPointStateType) :: point
    type(MaterialPointResponseType) :: response
    type(MaterialPointTangentType) :: tangent
    type(LoadStateType) :: load
    real(dp), allocatable :: gp(:),gw(:)
    real(dp) :: coords(3,20),shape(20),gradient(3,20),phi(20),eta_n(20),eta_p(20),theta(20)
    real(dp) :: weight,det,x,jn,jp,qje,qdrift,qdiff,peak,err(8),volume,exact_gradient
    integer :: element,i,j,k,node
    call independent_gauss_rule(5,gp,gw,status)
    if (.not. status%is_ok()) return
    call analytic_constants(slope,jn,jp,qje,qdrift,qdiff)
    peak=qje*LENGTH**2/(8.0_dp*KAPPA)
    errors=0.0_dp; point_errors=0.0_dp; volume=0.0_dp; drift_integral=0.0_dp
    load%voltage_factor=1.0_dp; load%joule_factor=1.0_dp
    do element=1,mesh%number_of_elements()
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      do node=1,20
        phi(node)=state%x(global_dof(mesh%connectivity(node,element),F_PHI))
        eta_n(node)=state%x(global_dof(mesh%connectivity(node,element),F_ETA_N))
        eta_p(node)=state%x(global_dof(mesh%connectivity(node,element),F_ETA_P))
        theta(node)=state%x(global_dof(mesh%connectivity(node,element),F_THETA))
      end do
      do k=1,5
        do j=1,5
          do i=1,5
            call hex20_geometry(coords,gp(i),gp(j),gp(k),shape,gradient,det,status)
            if (.not. status%is_ok()) return
            weight=gw(i)*gw(j)*gw(k)*det; volume=volume+weight
            x=dot_product(shape,coords(1,:))
            point=MaterialPointStateType()
            point%electric_field=-matmul(gradient,phi)
            point%eta_n=dot_product(shape,eta_n); point%eta_p=dot_product(shape,eta_p)
            point%eta_n_gradient=matmul(gradient,eta_n); point%eta_p_gradient=matmul(gradient,eta_p)
            point%temperature=T0+dot_product(shape,theta)
            point%temperature_gradient=matmul(gradient,theta)
            call evaluate_material_point(material,point,load,response,tangent,status)
            if (.not. status%is_ok()) return
            exact_gradient=qje*(LENGTH-2.0_dp*x)/(2.0_dp*KAPPA)
            err(1)=((dot_product(shape,phi)+VOLTAGE*x/LENGTH)/VOLTAGE)**2
            err(2)=sum((point%electric_field-[VOLTAGE/LENGTH,0.0_dp,0.0_dp])**2)*(LENGTH/VOLTAGE)**2
            err(3)=((point%eta_n-slope*x/LENGTH)/ETA_SCALE)**2
            err(4)=((point%eta_p-slope*x/LENGTH)/ETA_SCALE)**2
            err(5)=sum((point%eta_n_gradient-[slope/LENGTH,0.0_dp,0.0_dp])**2)*(LENGTH/ETA_SCALE)**2
            err(6)=sum((point%eta_p_gradient-[slope/LENGTH,0.0_dp,0.0_dp])**2)*(LENGTH/ETA_SCALE)**2
            err(7)=((dot_product(shape,theta)-qje*x*(LENGTH-x)/(2.0_dp*KAPPA))/peak)**2
            err(8)=sum((point%temperature_gradient-[exact_gradient,0.0_dp,0.0_dp])**2)*(LENGTH/peak)**2
            errors=errors+weight*err
            point_errors(1)=point_errors(1)+weight*sum((response%current_n-[jn,0.0_dp,0.0_dp])**2)/jn**2
            point_errors(2)=point_errors(2)+weight*sum((response%current_p-[jp,0.0_dp,0.0_dp])**2)/jp**2
            point_errors(3)=point_errors(3)+weight*((response%joule_heat-qje)/qje)**2
            ! LQ currents use reference densities, not response%n/response%p.
            drift_integral=drift_integral+weight*ELEMENTARY_CHARGE*DENSITY*(MU_N+MU_P)* &
              sum(point%electric_field**2)
          end do
        end do
      end do
    end do
    errors=sqrt(errors/volume); point_errors=sqrt(point_errors/volume)
  end subroutine field_errors

  subroutine verify_nt_material(passed,status)
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(MaterialType) :: material,old_material
    type(MaterialPointStateType) :: point
    type(MaterialPointResponseType) :: response,old_response
    type(MaterialPointTangentType) :: tangent,old_tangent
    type(LoadStateType) :: load
    real(dp) :: mn,mp,dn,dpv,n,p,vthermal,e(3),gn(3),gp(3),jn(3),jp(3)
    real(dp) :: qn,qp,qdrift,qexpected,eq_e(3),eq_gn(3),eq_gp(3),eq_t,errors(8),local_error
    integer :: mode
    call define_material(material)
    material%linear_carriers=.false.; material%frozen_transport=.false.
    material%beta_n=1.35_dp; material%beta_p=1.10_dp
    old_material=material; old_material%heat_source_current_work=.false.
    load%voltage_factor=1.0_dp; load%joule_factor=1.0_dp
    point=MaterialPointStateType()
    point%temperature=347.0_dp; point%eta_n=0.13_dp; point%eta_p=-0.07_dp
    e=[1.3e4_dp,-0.7e4_dp,0.4e4_dp]
    point%electric_field=e
    mn=MU_N*(point%temperature/T0)**(-material%beta_n)
    mp=MU_P*(point%temperature/T0)**(-material%beta_p)
    vthermal=BOLTZMANN_CONSTANT*point%temperature/ELEMENTARY_CHARGE
    dn=mn*vthermal; dpv=mp*vthermal
    n=DENSITY*exp(point%eta_n); p=DENSITY*exp(point%eta_p)
    qdrift=ELEMENTARY_CHARGE*(n*mn+p*mp)*sum(e**2)
    passed=.true.
    do mode=1,2
      if (mode == 1) then
        gn=-e/vthermal; gp=e/vthermal
      else
        gn=[-1.1e5_dp,0.9e5_dp,0.6e5_dp]
        gp=[0.8e5_dp,-1.4e5_dp,0.5e5_dp]
      end if
      point%eta_n_gradient=gn; point%eta_p_gradient=gp
      call evaluate_material_point(material,point,load,response,tangent,status)
      if (.not. status%is_ok()) return
      call evaluate_material_point(old_material,point,load,old_response,old_tangent,status)
      if (.not. status%is_ok()) return
      ! Independent scalar-isotropic closed formulas; no production transport outputs reused.
      jn=ELEMENTARY_CHARGE*n*(mn*e+dn*gn)
      jp=ELEMENTARY_CHARGE*p*(mp*e-dpv*gp)
      qn=ELEMENTARY_CHARGE*n*(mn*sum(e**2)+dn*dot_product(e,gn))
      qp=ELEMENTARY_CHARGE*p*(mp*sum(e**2)-dpv*dot_product(e,gp))
      qexpected=qn+qp
      eq_e=ELEMENTARY_CHARGE*(2.0_dp*(n*mn+p*mp)*e+n*dn*gn-p*dpv*gp)
      eq_gn=ELEMENTARY_CHARGE*n*dn*e; eq_gp=-ELEMENTARY_CHARGE*p*dpv*e
      eq_t=ELEMENTARY_CHARGE/point%temperature*( &
        -n*material%beta_n*mn*sum(e**2)+n*(1.0_dp-material%beta_n)*dn*dot_product(e,gn) &
        -p*material%beta_p*mp*sum(e**2)-p*(1.0_dp-material%beta_p)*dpv*dot_product(e,gp))
      errors(1)=maxval(abs(response%current_n-jn))/(ELEMENTARY_CHARGE*n*mn*sqrt(sum(e**2)))
      errors(2)=maxval(abs(response%current_p-jp))/(ELEMENTARY_CHARGE*p*mp*sqrt(sum(e**2)))
      errors(3)=abs(response%joule_heat-qexpected)/qdrift
      errors(4)=maxval(abs(tangent%joule_E-eq_e))/(qdrift/sqrt(sum(e**2)))
      errors(5)=max(abs(tangent%joule_eta_n-qn),abs(tangent%joule_eta_p-qp))/qdrift
      errors(6)=maxval(abs(tangent%joule_grad_eta_n-eq_gn))/maxval(abs(eq_gn))
      errors(7)=maxval(abs(tangent%joule_grad_eta_p-eq_gp))/maxval(abs(eq_gp))
      errors(8)=abs(tangent%joule_T-eq_t)/(qdrift/point%temperature)
      local_error=maxval(errors)
      passed=passed .and. local_error < FIELD_TOL .and. &
        abs(old_response%joule_heat/qdrift-1.0_dp) < POWER_TOL
      if (mode == 1) then
        local_error=max(abs(response%joule_heat)/qdrift, &
          maxval(abs(response%current_n))/(ELEMENTARY_CHARGE*n*mn*sqrt(sum(e**2))), &
          maxval(abs(response%current_p))/(ELEMENTARY_CHARGE*p*mp*sqrt(sum(e**2))))
        passed=passed .and. local_error < FIELD_TOL .and. old_response%joule_heat > 0.0_dp
        call metric('NT_local_equilibrium','cancellation_scaled_error',local_error)
        call metric('NT_local_equilibrium','old_drift_Q_W_m3',old_response%joule_heat)
        call metric('NT_local_equilibrium','JE_Q_W_m3',response%joule_heat)
        call metric('NT_local_equilibrium','max_independent_tangent_error',maxval(errors(4:8)))
      else
        call metric('NT_nonzero_gradient','max_independent_response_error',maxval(errors(1:3)))
        call metric('NT_nonzero_gradient','max_independent_tangent_error',maxval(errors(4:8)))
      end if
    end do
    write(*,'(a)') 'jegrad_NT_scope=local_material_only_not_global_Poisson_equilibrium_or_thermoelectric_closure'
    write(*,'(a,l1)') 'jegrad_NT_local_checks_passed=',passed
  end subroutine verify_nt_material

  subroutine metric(label,name,value)
    character(len=*), intent(in) :: label,name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') 'jegrad_'//trim(label)//'_'//trim(name)//'=',value
  end subroutine metric
end module je_gradient_check_m

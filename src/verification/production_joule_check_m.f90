module production_joule_check_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use field_layout_m, only : NDOF_NODE, F_U1, F_U2, F_U3, F_PHI, F_ETA_N, &
    F_ETA_P, F_THETA, global_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use model_types_m, only : MaterialType, StateType, LoadStateType, &
    NonlinearOptionsType, ContinuationOptionsType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use hex20_shape_m, only : hex20_natural_coordinates, hex20_geometry
  use mesh_checks_m, only : validate_mesh
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use newton_solver_m, only : NewtonReportType, solve_single_load_newton
  use continuation_solver_m, only : ContinuationReportType, solve_load_continuation
  use independent_flux_audit_m, only : IndependentFluxReportType, &
    integrate_independent_fluxes, independent_gauss_rule
  implicit none
  private

  real(dp), parameter :: LENGTH=20.0e-6_dp, WIDTH=20.0e-6_dp, THICKNESS=1.0e-6_dp
  real(dp), parameter :: FIELD_SCALES(7)=[1.0e-10_dp,1.0e-10_dp,1.0e-10_dp, &
    3.0_dp,0.1_dp,0.1_dp,10.0_dp]
  real(dp), parameter :: FIELD_TOL=1.0e-8_dp, FLUX_TOL=1.0e-8_dp
  integer, parameter :: COUNTS(4)=[1,2,4,8]
  public :: verify_production_joule_a1j

contains

  subroutine verify_production_joule_a1j(status,current_work)
    type(SolverStatusType), intent(out) :: status
    logical, optional, intent(in) :: current_work
    integer :: profile,level,total
    logical :: passed,one_pass,use_current_work
    character(len=48) :: label

    call clear_status(status)
    use_current_work=.false.
    if (present(current_work)) use_current_work=current_work
    passed=.true.; total=0
    if (use_current_work) then
      write(*,'(a)') 'verification_id=A1JJE_production_joule_independent_flux'
      write(*,'(a)') 'a1j_heat_source=current_work'
    else
      write(*,'(a)') 'verification_id=A1J_production_joule_independent_flux'
    end if
    write(*,'(a)') 'a1j_scope=production_material_assembly_newton_continuation_not_native_run_profile'
    write(*,'(a)') 'a1j_interpolation=linear_phi_quadratic_theta_exactly_representable'
    do profile=1,2
      do level=1,size(COUNTS)
        write(label,'(a,i0,a,i0)') 'p',profile,'_nx',COUNTS(level)
        call solve_case(trim(label),profile,COUNTS(level),1,1,3.0_dp,1.0_dp, &
          level == 2,one_pass,status,use_current_work)
        if (.not. status%is_ok()) return
        passed=passed .and. one_pass; total=total+1
      end do
    end do
    call solve_case('reverse_voltage',1,2,1,1,-3.0_dp,1.0_dp,.false.,one_pass,status,use_current_work)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass; total=total+1
    call solve_case('zero_voltage',1,2,1,1,0.0_dp,1.0_dp,.false.,one_pass,status,use_current_work)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass; total=total+1
    call solve_case('half_voltage',1,2,1,1,1.5_dp,1.0_dp,.false.,one_pass,status,use_current_work)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass; total=total+1
    call solve_case('joule_off',1,2,1,1,3.0_dp,0.0_dp,.false.,one_pass,status,use_current_work)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass; total=total+1
    call solve_case('transverse_refinement',1,2,2,2,3.0_dp,1.0_dp,.false.,one_pass,status,use_current_work)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass; total=total+1
    call verify_flux_patch(one_pass,status,use_current_work)
    if (.not. status%is_ok()) return
    passed=passed .and. one_pass
    write(*,'(a,i0)') 'a1j_solved_cases=',total
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'A1J exceeded a frozen field/flux/negative-control criterion.')
    end if
  end subroutine verify_production_joule_a1j

  subroutine solve_case(label,profile,nx,ny,nz,voltage,joule_factor,extra_checks,passed,status,current_work)
    character(len=*), intent(in) :: label
    integer, intent(in) :: profile,nx,ny,nz
    real(dp), intent(in) :: voltage,joule_factor
    logical, intent(in) :: extra_checks,current_work
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state,second_state,bad_state
    type(LoadStateType) :: load,start_load
    type(NonlinearOptionsType) :: options
    type(ContinuationOptionsType) :: continuation_options
    type(NewtonReportType) :: newton,second_newton
    type(ContinuationReportType) :: continuation
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: boundaries(:)
    type(IndependentFluxReportType) :: flux,flux4,bad_flux
    real(dp) :: sigma_n,sigma_p,power,power_scale,current_n,current_p,temperature_peak
    real(dp) :: errors(6),flux_errors(9),quad_error,geometry_error,state_difference
    real(dp) :: minimum_det,volume,xi,bad_heat_defect,bad_current_defect
    integer :: node,expected_exterior,expected_internal

    call clear_status(status)
    passed=.false.
    call build_box(nx,ny,nz,mesh)
    call define_material(profile,materials(1))
    materials(1)%heat_source_current_work=current_work
    call validate_mesh(mesh,minimum_det,volume,status)
    if (.not. status%is_ok()) return
    call define_boundaries(mesh,voltage,boundaries)
    call initialize_sparse_context(mesh,context,status)
    if (.not. status%is_ok()) return
    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()))
    state%x=0.0_dp
    load%voltage_factor=1.0_dp; load%joule_factor=joule_factor
    options%max_newton=30
    options%residual_tolerance=1.0e-11_dp
    options%update_tolerance=1.0e-11_dp
    call solve_single_load_newton(mesh,materials,load,context,boundaries,FIELD_SCALES, &
      options,state,newton,status)
    if (.not. status%is_ok()) then
      write(*,'(a)') 'a1j_failed_case='//trim(label)
      call metric(label,'failed_residual',newton%final_residual_rms)
      return
    end if
    call field_errors(mesh,materials(1),state,voltage,joule_factor,errors,status)
    if (.not. status%is_ok()) return
    call integrate_independent_fluxes(mesh,materials,state,load,5,flux,status)
    if (.not. status%is_ok()) return
    call integrate_independent_fluxes(mesh,materials,state,load,4,flux4,status)
    if (.not. status%is_ok()) return
    sigma_n=ELEMENTARY_CHARGE*materials(1)%reference_n*materials(1)%mobility_n0(1,1)
    sigma_p=ELEMENTARY_CHARGE*materials(1)%reference_p*materials(1)%mobility_p0(1,1)
    current_n=sigma_n*voltage/LENGTH*WIDTH*THICKNESS
    current_p=sigma_p*voltage/LENGTH*WIDTH*THICKNESS
    power=(sigma_n+sigma_p)*voltage**2/LENGTH*WIDTH*THICKNESS
    power_scale=(sigma_n+sigma_p)*max(voltage**2,1.0_dp)/LENGTH*WIDTH*THICKNESS
    temperature_peak=joule_factor*(sigma_n+sigma_p)*voltage**2/ &
      (8.0_dp*materials(1)%thermal_conductivity(1,1))
    flux_errors(1)=maxval(abs(flux%electron_out-[-current_n,current_n,0.0_dp]))/ &
      (sigma_n*max(abs(voltage),1.0_dp)*WIDTH*THICKNESS/LENGTH)
    flux_errors(2)=maxval(abs(flux%hole_out-[-current_p,current_p,0.0_dp]))/ &
      (sigma_p*max(abs(voltage),1.0_dp)*WIDTH*THICKNESS/LENGTH)
    flux_errors(3)=maxval(abs(flux%heat_out-[0.5_dp*joule_factor*power, &
      0.5_dp*joule_factor*power,0.0_dp]))/power_scale
    flux_errors(4)=abs(flux%joule-joule_factor*power)/power_scale
    flux_errors(5)=abs(sum(flux%heat_out)-flux%joule)/power_scale
    flux_errors(6)=abs(flux%terminal_power-power)/power_scale
    flux_errors(7)=abs(flux%terminal_power-flux%total_current_work)/power_scale
    flux_errors(8)=abs(flux%electric_displacement_out-flux%charge)/ &
      max(materials(1)%permittivity(1,1)*max(abs(voltage),1.0_dp)/LENGTH*WIDTH*THICKNESS, &
          flux%absolute_displacement_out+flux%absolute_charge)
    flux_errors(9)=flux%absolute_diffusion_work/power_scale
    quad_error=max(abs(flux%volume-flux4%volume)/(LENGTH*WIDTH*THICKNESS), &
      abs(flux%joule-flux4%joule)/power_scale, &
      maxval(abs(flux%heat_out-flux4%heat_out))/power_scale)
    geometry_error=max(abs(flux%volume/(LENGTH*WIDTH*THICKNESS)-1.0_dp), &
      abs(flux%area/(2.0_dp*(LENGTH*WIDTH+WIDTH*THICKNESS+LENGTH*THICKNESS))-1.0_dp), &
      sqrt(sum(flux%area_normal**2))/flux%area)
    expected_exterior=2*(nx*ny+ny*nz+nz*nx)
    expected_internal=(nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1)
    passed=newton%converged .and. newton%maximum_perturbed_pivots == 0 .and. &
      maxval(errors) < FIELD_TOL .and. maxval(flux_errors) < FLUX_TOL .and. &
      quad_error < 1.0e-10_dp .and. geometry_error < 1.0e-10_dp .and. &
      flux%exterior_faces == expected_exterior .and. flux%internal_faces == expected_internal
    write(*,'(a,i0)') 'a1j_'//trim(label)//'_elements=',mesh%number_of_elements()
    write(*,'(a,i0)') 'a1j_'//trim(label)//'_newton_iterations=',newton%iterations
    call metric(label,'phi_L2_scaled',errors(1))
    call metric(label,'E_L2_scaled',errors(2))
    call metric(label,'theta_L2_scaled',errors(3))
    call metric(label,'grad_theta_L2_scaled',errors(4))
    call metric(label,'max_relative_n_error',errors(5))
    call metric(label,'max_relative_p_error',errors(6))
    call metric(label,'max_flux_error',maxval(flux_errors))
    call metric(label,'heat_balance_defect',flux_errors(5))
    call metric(label,'quadrature_4_vs_5',quad_error)
    call metric(label,'geometry_error',geometry_error)
    call metric(label,'electron_right_A',flux%electron_out(2))
    call metric(label,'hole_right_A',flux%hole_out(2))
    call metric(label,'joule_W',flux%joule)
    call metric(label,'heat_left_W',flux%heat_out(1))
    call metric(label,'heat_right_W',flux%heat_out(2))
    call metric(label,'terminal_power_W',flux%terminal_power)
    call metric(label,'analytic_peak_K',temperature_peak)
    call metric(label,'numerical_peak_K',maxval(state%x(F_THETA::NDOF_NODE)))

    if (extra_checks) then
      allocate(second_state%x(size(state%x)))
      do node=1,mesh%number_of_nodes()
        second_state%x(global_dof(node,F_U1):global_dof(node,F_THETA))=0.0_dp
        second_state%x(global_dof(node,F_PHI))=0.8_dp*state%x(global_dof(node,F_PHI))
        second_state%x(global_dof(node,F_THETA))=0.4_dp*state%x(global_dof(node,F_THETA))
      end do
      call solve_single_load_newton(mesh,materials,load,context,boundaries,FIELD_SCALES, &
        options,second_state,second_newton,status)
      if (.not. status%is_ok()) return
      state_difference=scaled_difference(state,second_state)
      passed=passed .and. state_difference < FIELD_TOL .and. second_newton%converged
      call metric(label,'dual_start_difference',state_difference)
      second_state%x=0.0_dp
      start_load%voltage_factor=0.0_dp; start_load%joule_factor=0.0_dp
      continuation_options%maximum_newton_per_attempt=30
      call solve_load_continuation(mesh,materials,context,boundaries,FIELD_SCALES,options, &
        continuation_options,load,second_state,start_load,continuation,status)
      if (.not. status%is_ok()) return
      state_difference=scaled_difference(state,second_state)
      passed=passed .and. state_difference < FIELD_TOL .and. continuation%converged
      call metric(label,'continuation_difference',state_difference)
      allocate(bad_state%x(size(state%x)))
      bad_state%x=state%x
      bad_state%x(F_THETA::NDOF_NODE)=1.1_dp*bad_state%x(F_THETA::NDOF_NODE)
      call integrate_independent_fluxes(mesh,materials,bad_state,load,5,bad_flux,status)
      if (.not. status%is_ok()) return
      bad_heat_defect=abs(sum(bad_flux%heat_out)-bad_flux%joule)/power_scale
      bad_state%x=state%x
      do node=1,mesh%number_of_nodes()
        xi=mesh%coordinates(1,node)/LENGTH
        bad_state%x(global_dof(node,F_PHI))=bad_state%x(global_dof(node,F_PHI))+ &
          0.2_dp*voltage*xi*(1.0_dp-xi)
      end do
      call integrate_independent_fluxes(mesh,materials,bad_state,load,5,bad_flux,status)
      if (.not. status%is_ok()) return
      bad_current_defect=abs(sum(bad_flux%electron_out))/max(abs(current_n),tiny(1.0_dp))
      passed=passed .and. bad_heat_defect > 0.05_dp .and. bad_current_defect > 0.1_dp
      call metric(label,'negative_heat_control_defect',bad_heat_defect)
      call metric(label,'negative_current_control_defect',bad_current_defect)
    end if
    write(*,'(a,l1)') 'a1j_'//trim(label)//'_passed=',passed
  end subroutine solve_case

  subroutine define_material(profile,material)
    integer, intent(in) :: profile
    type(MaterialType), intent(out) :: material
    integer :: i
    material=MaterialType()
    material%name='A1J_SYNTHETIC_TWO_CARRIER'
    material%reference_temperature=300.0_dp
    material%reference_n=1.0e21_dp; material%reference_p=8.0e20_dp
    do i=1,6
      material%elastic(i,i)=1.0e11_dp
    end do
    do i=1,3
      material%permittivity(i,i)=8.0e-11_dp
      material%thermal_conductivity(i,i)=10.0_dp
      material%mobility_n0(i,i)=0.05_dp
      material%mobility_p0(i,i)=0.02_dp
    end do
    if (profile == 2) then
      ! Only transport magnitudes and minority scale are borrowed from the candidate.
      ! No actual ZnO piezoelectric, pyroelectric or thermal-expansion claim.
      material%name='A1J_ELECTRON_DOMINANT_SCALE_TEST'
      material%reference_n=6.0e22_dp
      material%reference_p=(2.9e-4_dp)**2/material%reference_n
      material%mobility_n0=0.0_dp; material%mobility_p0=0.0_dp
      material%thermal_conductivity=0.0_dp
      do i=1,3
        material%mobility_n0(i,i)=0.0205_dp
        material%mobility_p0(i,i)=0.005_dp
        material%thermal_conductivity(i,i)=60.0_dp
      end do
    end if
    material%donor_concentration=material%reference_n-material%reference_p
    material%acceptor_concentration=0.0_dp
    material%intrinsic_concentration=sqrt(material%reference_n*material%reference_p)
    ! Defaults leave e,b,pyro and beta exactly zero; D(T) still follows Einstein.
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

  subroutine define_boundaries(mesh,voltage,boundaries)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: voltage
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
    call set_boundary(boundaries(5),F_PHI,right,-voltage)
    call set_boundary(boundaries(6),F_ETA_N,left,0.0_dp)
    call set_boundary(boundaries(7),F_ETA_N,right,0.0_dp)
    call set_boundary(boundaries(8),F_ETA_P,left,0.0_dp)
    call set_boundary(boundaries(9),F_ETA_P,right,0.0_dp)
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

  subroutine field_errors(mesh,material,state,voltage,joule_factor,errors,status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    real(dp), intent(in) :: voltage,joule_factor
    real(dp), intent(out) :: errors(6)
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: gp(:),gw(:)
    real(dp) :: coords(3,20),shape(20),gradient(3,20),phi(20),theta(20)
    real(dp) :: weight,det,x,exact_phi,exact_theta,exact_gradient,coefficient
    real(dp) :: sigma,phi_scale,theta_scale,grad_value(3),err(4),volume
    integer :: element,i,j,k,node
    call independent_gauss_rule(5,gp,gw,status)
    errors=0.0_dp; volume=0.0_dp
    sigma=ELEMENTARY_CHARGE*(material%reference_n*material%mobility_n0(1,1)+ &
      material%reference_p*material%mobility_p0(1,1))
    coefficient=joule_factor*sigma*voltage**2/(2.0_dp*material%thermal_conductivity(1,1)*LENGTH**2)
    phi_scale=max(abs(voltage),1.0_dp)
    theta_scale=max(coefficient*LENGTH**2/4.0_dp,1.0_dp)
    do element=1,mesh%number_of_elements()
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      do node=1,20
        phi(node)=state%x(global_dof(mesh%connectivity(node,element),F_PHI))
        theta(node)=state%x(global_dof(mesh%connectivity(node,element),F_THETA))
      end do
      do k=1,5
        do j=1,5
          do i=1,5
            call hex20_geometry(coords,gp(i),gp(j),gp(k),shape,gradient,det,status)
            if (.not. status%is_ok()) return
            weight=gw(i)*gw(j)*gw(k)*det; volume=volume+weight
            x=dot_product(shape,coords(1,:))
            exact_phi=-voltage*x/LENGTH
            exact_theta=coefficient*x*(LENGTH-x)
            exact_gradient=coefficient*(LENGTH-2.0_dp*x)
            grad_value=matmul(gradient,phi)-[-voltage/LENGTH,0.0_dp,0.0_dp]
            err(1)=((dot_product(shape,phi)-exact_phi)/phi_scale)**2
            err(2)=sum(grad_value**2)*(LENGTH/phi_scale)**2
            err(3)=((dot_product(shape,theta)-exact_theta)/theta_scale)**2
            grad_value=matmul(gradient,theta)-[exact_gradient,0.0_dp,0.0_dp]
            err(4)=sum(grad_value**2)*(LENGTH/theta_scale)**2
            errors(1:4)=errors(1:4)+weight*err
          end do
        end do
      end do
    end do
    errors(1:4)=sqrt(errors(1:4)/volume)
    errors(5)=maxval(abs(exp(state%x(F_ETA_N::NDOF_NODE))-1.0_dp))
    errors(6)=maxval(abs(exp(state%x(F_ETA_P::NDOF_NODE))-1.0_dp))
  end subroutine field_errors

  subroutine verify_flux_patch(passed,status,current_work)
    logical, intent(out) :: passed
    logical, intent(in) :: current_work
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(IndependentFluxReportType) :: flux
    real(dp) :: x(3),expected_heat,expected_joule,error_heat,error_joule,efield(3)
    integer :: node
    call build_box(2,2,2,mesh)
    call define_material(1,materials(1))
    materials(1)%heat_source_current_work=current_work
    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()))
    state%x=0.0_dp
    efield=[1.0e3_dp,2.0e3_dp,3.0e3_dp]
    do node=1,mesh%number_of_nodes()
      x=mesh%coordinates(:,node)
      state%x(global_dof(node,F_PHI))=-dot_product(efield,x)
      state%x(global_dof(node,F_THETA))=sum((x/[LENGTH,WIDTH,THICKNESS])**2)
    end do
    load%voltage_factor=1.0_dp; load%joule_factor=1.0_dp
    call integrate_independent_fluxes(mesh,materials,state,load,5,flux,status)
    if (.not. status%is_ok()) return
    expected_heat=-2.0_dp*10.0_dp*sum(1.0_dp/[LENGTH,WIDTH,THICKNESS]**2)*LENGTH*WIDTH*THICKNESS
    expected_joule=ELEMENTARY_CHARGE*(1.0e21_dp*0.05_dp+8.0e20_dp*0.02_dp)* &
      sum(efield**2)*LENGTH*WIDTH*THICKNESS
    error_heat=abs(sum(flux%heat_out)/expected_heat-1.0_dp)
    error_joule=abs(flux%joule/expected_joule-1.0_dp)
    passed=error_heat < 1.0e-10_dp .and. error_joule < 1.0e-10_dp .and. &
      flux%exterior_faces == 24 .and. flux%internal_faces == 12
    call metric('flux_patch','three_axis_heat_flux_error',error_heat)
    call metric('flux_patch','three_axis_joule_error',error_joule)
    write(*,'(a,l1)') 'a1j_flux_patch_passed=',passed
  end subroutine verify_flux_patch

  real(dp) function scaled_difference(a,b) result(value)
    type(StateType), intent(in) :: a,b
    integer :: field
    value=0.0_dp
    do field=1,NDOF_NODE
      value=max(value,maxval(abs(a%x(field::NDOF_NODE)-b%x(field::NDOF_NODE)))/FIELD_SCALES(field))
    end do
  end function scaled_difference

  subroutine metric(label,name,value)
    character(len=*), intent(in) :: label,name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') 'a1j_'//trim(label)//'_'//trim(name)//'=',value
  end subroutine metric
end module production_joule_check_m

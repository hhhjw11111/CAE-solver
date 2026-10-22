module f1_case_runner_m
  use, intrinsic :: iso_fortran_env, only : output_unit
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use field_layout_m, only : NDOF_NODE,F_PHI,F_ETA_N,F_ETA_P,F_THETA,global_dof
  use input_schema_m, only : CaseConfigType
  use model_types_m, only : MaterialType,StateType,LoadStateType,DiagnosticsType
  use mesh_types_m, only : MeshType
  use native_input_m, only : read_native_case
  use mesh_checks_m, only : validate_mesh
  use zno_candidate_m, only : define_zno_candidate
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_INPUT_ERROR,STATUS_VERIFICATION_FAILED
  use dirichlet_m, only : DirichletBCType
  use sparse_pattern_m, only : SparseAssemblyContextType,initialize_sparse_context
  use newton_solver_m, only : NewtonReportType,solve_single_load_newton
  use continuation_solver_m, only : ContinuationReportType,solve_load_continuation
  use balance_audit_m, only : audit_global_balances
  use independent_flux_audit_m, only : IndependentFluxReportType,integrate_independent_fluxes,independent_gauss_rule
  use field_metrics_m, only : HoleFieldMetricsType,evaluate_hole_field_metrics,weighted_percentile
  use hex20_shape_m, only : hex20_geometry
  use element_kernel_m, only : build_operators,unpack_point_state
  use material_point_m, only : MaterialPointStateType,MaterialPointResponseType,MaterialPointTangentType, &
    evaluate_material_point,MATERIAL_INPUT_SIZE,MATERIAL_RESPONSE_SIZE
  use dof_map_m, only : HEX20_NDOF,build_element_dof_map
  implicit none
  private
  public :: run_f1_case
contains
  subroutine run_f1_case(case_path,directory,status)
    character(len=*), intent(in) :: case_path,directory
    type(SolverStatusType), intent(out) :: status
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state
    type(LoadStateType) :: load,target
    type(DirichletBCType), allocatable :: bc(:),stage_bc(:)
    type(SparseAssemblyContextType) :: context
    type(NewtonReportType) :: newton
    type(ContinuationReportType) :: continuation
    type(DiagnosticsType) :: diagnostics
    type(IndependentFluxReportType) :: flux,flux4
    type(HoleFieldMetricsType) :: ring
    character(len=2) :: model
    real(dp) :: minimum_det,volume,balance(4),bc_error,max_balance,max_residual,max_update
    real(dp) :: heat_error,jn_error,jp_error,quad_error,power_scale,current_scale,extra(8)
    integer :: unit_id,io,stage,point,total_iterations
    logical :: accepted
    call clear_status(status)
    call read_native_case(case_path,config,mesh,status)
    if (.not. status%is_ok()) return
    select case(trim(config%material_profile))
    case('ZNO_F1_LQ');model='LQ'
    case('ZNO_F1_N0');model='N0'
    case('ZNO_F1_NT');model='NT'
    case default
      call set_status(status,STATUS_INPUT_ERROR,'run-f1 requires an explicit ZNO_F1 candidate profile.');return
    end select
    open(newunit=unit_id,file=trim(directory)//'/metrics.txt',status='new',action='write',iostat=io)
    if (io/=0) then
      call set_status(status,STATUS_INPUT_ERROR,'F1 requires an existing NEW result directory.');return
    end if
    write(unit_id,'(a)') 'case_id='//trim(config%case_id)
    write(unit_id,'(a)') 'model='//model
    write(unit_id,'(a)') 'material_evidence=NUMERICAL_CANDIDATE_NOT_EXPERIMENTALLY_CALIBRATED'
    write(unit_id,'(a)') 'source=ZNO_300K_CANDIDATE_V2; plan=F1_20260906'
    write(unit_id,'(a)') 'heat_source=DRIFT_CONDUCTIVITY_NOT_TOTAL_CURRENT_WORK'
    flush(unit_id)
    call validate_mesh(mesh,minimum_det,volume,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call define_zno_candidate(model,materials(1))
    config%continuation_options%heat_balance_reference_power=ELEMENTARY_CHARGE*( &
      materials(1)%reference_n*materials(1)%mobility_n0(1,1)+ &
      materials(1)%reference_p*materials(1)%mobility_p0(1,1))*config%target_voltage**2*1.0e-6_dp
    write(unit_id,'(a,es24.16e3)') 'zero_source_heat_reference_W=', &
      config%continuation_options%heat_balance_reference_power
    call define_boundaries(mesh,config%target_voltage,bc,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call initialize_sparse_context(mesh,context,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    write(output_unit,'(a,i0,a,i0)') 'f1_nodes=',mesh%number_of_nodes(),', elements=',mesh%number_of_elements()
    flush(output_unit)
    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()));state%x=0.0_dp
    config%nonlinear_options%verbose=.true.
    max_balance=0.0_dp;max_residual=0.0_dp;max_update=0.0_dp;total_iterations=0
    if (model=='LQ') then
      write(unit_id,'(a)') 'solution_mode=THREE_STAGE_LQ_NO_RETURN_FEEDBACK'
      do stage=1,3
        call stage_boundaries(mesh,bc,state,stage,stage_bc)
        load%voltage_factor=1.0_dp;load%joule_factor=1.0_dp
        if (stage==1) load%joule_factor=0.0_dp
        write(output_unit,'(a,i0)') 'lq_stage=',stage
        call solve_single_load_newton(mesh,materials,load,context,stage_bc,config%field_scales, &
          config%nonlinear_options,state,newton,status)
        if (.not. status%is_ok()) then
          close(unit_id);return
        end if
        call audit_global_balances(mesh,materials,state,load,context,stage_bc,balance,bc_error,diagnostics,status)
        if (.not. status%is_ok()) then
          close(unit_id);return
        end if
        max_balance=max(max_balance,maxval(balance));max_residual=max(max_residual,newton%final_residual_rms)
        max_update=max(max_update,newton%final_update_rms);total_iterations=total_iterations+newton%iterations
        write(unit_id,'(a,i0,a,3(es24.16e3,a))') 'stage=',stage,',residual=',newton%final_residual_rms, &
          ',update=',newton%final_update_rms,',balance=',maxval(balance),','
        flush(unit_id)
      end do
    else
      load%voltage_factor=0.0_dp;load%joule_factor=0.0_dp
      target%voltage_factor=1.0_dp;target%joule_factor=1.0_dp
      call solve_load_continuation(mesh,materials,context,bc,config%field_scales,config%nonlinear_options, &
        config%continuation_options,target,state,load,continuation,status)
      if (.not. status%is_ok()) then
        close(unit_id);return
      end if
      max_balance=continuation%maximum_balance_defect
      max_residual=maxval(continuation%residual_rms(1:continuation%accepted_points))
      max_update=maxval(continuation%update_rms(1:continuation%accepted_points))
      total_iterations=continuation%total_newton_iterations
      write(unit_id,'(a)') 'continuation_columns=phase,voltage,joule,iterations,residual,update,max_balance'
      do point=1,continuation%accepted_points
        write(unit_id,'(i0,a,2(es24.16e3,a),i0,a,3(es24.16e3,a))') continuation%phase(point),',', &
          continuation%voltage_factor(point),',',continuation%joule_factor(point),',', &
          continuation%newton_iterations(point),',',continuation%residual_rms(point),',', &
          continuation%update_rms(point),',',maxval(continuation%balance(:,point)),','
        write(unit_id,'(a,i0,a,4(es24.16e3,","),2(es24.16e3,:,","))') 'balance_detail_point=',point,',', &
          continuation%balance(:,point),continuation%heat_free_sum(point),continuation%heat_absolute_sum(point)
      end do
    end if
    write(output_unit,'(a)') 'f1_solve=COMPLETE; independent_flux_and_fields=START'
    flush(output_unit)
    ! Save converged fields before postprocessing acceptance; no result is lost on a failed gate.
    call write_nodes(directory,mesh,materials(1),state,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call integrate_independent_fluxes(mesh,materials,state,load,5,flux,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call integrate_independent_fluxes(mesh,materials,state,load,4,flux4,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call evaluate_hole_field_metrics(mesh,materials,state,load,3.0e-6_dp,0.5e-6_dp,ring,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call extra_fields(directory,mesh,materials(1),state,load,extra,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    power_scale=max(abs(flux%joule),abs(sum(flux%heat_out)),1.7736095e-3_dp*1.0e-12_dp)
    current_scale=max(sum(abs(flux%electron_out(1:2))),tiny(1.0_dp))
    heat_error=abs(sum(flux%heat_out)-flux%joule)/power_scale
    jn_error=abs(sum(flux%electron_out))/current_scale
    jp_error=abs(sum(flux%hole_out))/max(sum(abs(flux%hole_out(1:2))),tiny(1.0_dp))
    quad_error=abs(flux%joule-flux4%joule)/max(abs(flux%joule),1.7736095e-15_dp)
    call metric('voltage_V',abs(config%target_voltage))
    call metric('minimum_jacobian',minimum_det)
    call metric('total_volume_m3',flux%volume)
    call metric('maximum_stage_balance',max_balance)
    call metric('maximum_accepted_residual',max_residual)
    call metric('maximum_accepted_update',max_update)
    call metric('joule_W',flux%joule)
    call metric('terminal_power_W',flux%terminal_power)
    call metric('heat_out_W',sum(flux%heat_out))
    call metric('heat_balance_relative',heat_error)
    call metric('electron_balance_relative',jn_error)
    call metric('hole_balance_diagnostic',jp_error)
    call metric('joule_quadrature_4_5_relative',quad_error)
    call metric('electron_right_A',flux%electron_out(2))
    call metric('hole_right_A',flux%hole_out(2))
    call metric('ring_volume_m3',ring%annulus_volume)
    call metric('ring_Q_mean_W_m3',ring%annulus_joule_mean)
    call metric('ring_Q_p95_W_m3',ring%annulus_joule_p95)
    call metric('ring_E_mean_V_m',ring%annulus_electric_field_mean)
    call metric('ring_abs_n_relative_mean',ring%annulus_n_relative_mean)
    call metric('ring_abs_p_relative_mean',ring%annulus_p_relative_mean)
    call metric('ring_theta_mean_K',ring%annulus_temperature_rise_mean)
    call metric('ring_displacement_mean_m',ring%annulus_displacement_mean)
    call metric('E_p95_V_m',extra(1));call metric('E_above_1e6_volume_fraction',extra(2))
    call metric('theta_max_K',extra(3));call metric('maximum_abs_n_relative',extra(4))
    call metric('maximum_abs_p_relative',extra(5));call metric('displacement_max_m',extra(6))
    call metric('minimum_mu_n_ratio',extra(7));call metric('stress_von_mises_max_Pa',extra(8))
    write(unit_id,'(a,i0)') 'total_newton_iterations=',total_iterations
    accepted=max_balance<1.0e-8_dp .and. max_residual<1.0e-9_dp .and. max_update<1.0e-9_dp .and. &
      heat_error<0.01_dp .and. jn_error<0.01_dp .and. quad_error<0.001_dp
    if (accepted) then
      write(unit_id,'(a)') 'case_acceptance=PASS_NOT_MESH_ACCEPTANCE'
      write(output_unit,'(a)') 'case_acceptance=PASS_NOT_MESH_ACCEPTANCE'
    else
      write(unit_id,'(a)') 'case_acceptance=FAIL'
      write(output_unit,'(a)') 'case_acceptance=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'F1 case failed a frozen solve/flux criterion; fields retained.')
    end if
    close(unit_id)
  contains
    subroutine metric(name,value)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      write(unit_id,'(a,es24.16e3)') trim(name)//'=',value
      write(output_unit,'(a,es24.16e3)') trim(name)//'=',value
    end subroutine metric
  end subroutine run_f1_case

  subroutine define_boundaries(mesh,voltage,bc,status)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: voltage
    type(DirichletBCType), allocatable, intent(out) :: bc(:)
    type(SolverStatusType), intent(out) :: status
    integer, allocatable :: nodes(:),left(:),right(:),a(:),b(:),c(:)
    real(dp), parameter :: tol=1.0e-15_dp
    integer :: i
    call clear_status(status)
    nodes=[(i,i=1,mesh%number_of_nodes())]
    left=pack(nodes,abs(mesh%coordinates(1,:)+10e-6_dp)<tol)
    right=pack(nodes,abs(mesh%coordinates(1,:)-10e-6_dp)<tol)
    a=pack(nodes,maxval(abs(mesh%coordinates-spread([-10e-6_dp,-10e-6_dp,-0.5e-6_dp],2,size(nodes))),dim=1)<tol)
    b=pack(nodes,maxval(abs(mesh%coordinates-spread([-10e-6_dp,10e-6_dp,-0.5e-6_dp],2,size(nodes))),dim=1)<tol)
    c=pack(nodes,maxval(abs(mesh%coordinates-spread([-10e-6_dp,-10e-6_dp,0.5e-6_dp],2,size(nodes))),dim=1)<tol)
    if (size(a)/=1 .or. size(b)/=1 .or. size(c)/=1 .or. size(left)<8 .or. size(right)<8) then
      call set_status(status,STATUS_INPUT_ERROR,'F1 mesh lacks prescribed contacts or rigid-body anchors.');return
    end if
    allocate(bc(11))
    bc(1)%field_id=1;bc(1)%node_ids=[a,b,c]
    bc(2)%field_id=2;bc(2)%node_ids=a
    bc(3)%field_id=3;bc(3)%node_ids=[a,b]
    do i=4,7
      bc(2*i-4)%field_id=i;bc(2*i-4)%node_ids=left
      bc(2*i-3)%field_id=i;bc(2*i-3)%node_ids=right
    end do
    bc(5)%value_end=voltage
  end subroutine define_boundaries

  subroutine stage_boundaries(mesh,bc,state,stage,stage_bc)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), intent(in) :: bc(:)
    type(StateType), intent(in) :: state
    integer, intent(in) :: stage
    type(DirichletBCType), allocatable, intent(out) :: stage_bc(:)
    integer :: node,field,count_bc,k
    logical :: frozen(7)
    frozen=.false.
    if (stage==1) frozen(F_THETA)=.true.
    if (stage==2) frozen(1:6)=.true.
    if (stage==3) frozen(4:7)=.true.
    count_bc=size(bc)+count(frozen)*mesh%number_of_nodes()
    allocate(stage_bc(count_bc));stage_bc(1:size(bc))=bc;k=size(bc)
    do field=1,7
      if (.not. frozen(field)) cycle
      do node=1,mesh%number_of_nodes()
        k=k+1;stage_bc(k)%field_id=field;stage_bc(k)%node_ids=[node]
        stage_bc(k)%value_start=state%x(global_dof(node,field))
        stage_bc(k)%value_end=stage_bc(k)%value_start
      end do
    end do
  end subroutine stage_boundaries

  subroutine write_nodes(directory,mesh,material,state,status)
    character(len=*), intent(in) :: directory
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    type(SolverStatusType), intent(out) :: status
    integer :: unit_id,io,node
    real(dp) :: nratio,pratio
    call clear_status(status)
    open(newunit=unit_id,file=trim(directory)//'/nodes.csv',status='new',action='write',iostat=io)
    if(io/=0) then
      call set_status(status,STATUS_INPUT_ERROR,'Cannot create F1 nodes.csv.');return
    end if
    write(unit_id,'(a)') 'node,x_m,y_m,z_m,ux_m,uy_m,uz_m,phi_V,eta_n,eta_p,theta_K,n_m3,p_m3'
    do node=1,mesh%number_of_nodes()
      if (material%linear_carriers) then
        nratio=1.0_dp+state%x(global_dof(node,F_ETA_N));pratio=1.0_dp+state%x(global_dof(node,F_ETA_P))
      else
        nratio=exp(state%x(global_dof(node,F_ETA_N)));pratio=exp(state%x(global_dof(node,F_ETA_P)))
      end if
      write(unit_id,'(i0,12(",",es24.16e3))') node,mesh%coordinates(:,node), &
        state%x(global_dof(node,1):global_dof(node,7)),material%reference_n*nratio,material%reference_p*pratio
    end do
    close(unit_id)
  end subroutine write_nodes

  subroutine extra_fields(directory,mesh,material,state,load,extra,status)
    character(len=*), intent(in) :: directory
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    real(dp), intent(out) :: extra(8)
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: gp(:),gw(:),efields(:),weights(:)
    real(dp) :: coords(3,20),shape(20),grad(3,20),det,op(21,140),testop(140,21),raw(21),ue(140),x(3)
    real(dp) :: emagnitude,w,vol,highvol,disp(3),vm
    type(MaterialPointStateType) :: ps
    type(MaterialPointResponseType) :: response
    type(MaterialPointTangentType) :: tangent
    integer :: element,i,j,k,node,dofs(140),n,unit_id,io
    call independent_gauss_rule(5,gp,gw,status)
    if (.not. status%is_ok()) return
    allocate(efields(125*mesh%number_of_elements()),weights(125*mesh%number_of_elements()))
    open(newunit=unit_id,file=trim(directory)//'/field_samples.csv',status='new',action='write',iostat=io)
    if(io/=0) then
      call set_status(status,STATUS_INPUT_ERROR,'Cannot create F1 field samples.');return
    end if
    write(unit_id,'(a)') 'x_m,y_m,z_m,weight_m3,Ex_V_m,Ey_V_m,Ez_V_m,Q_W_m3,theta_K,n_over_n0,p_over_p0,'// &
      'ux_m,uy_m,uz_m,mu_n_over_mu0,sxx_Pa,syy_Pa,szz_Pa,syz_Pa,szx_Pa,sxy_Pa'
    extra=0.0_dp;extra(7)=huge(1.0_dp);n=0;vol=0.0_dp;highvol=0.0_dp
    do element=1,mesh%number_of_elements()
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      call build_element_dof_map(mesh%connectivity(:,element),dofs,status)
      if (.not. status%is_ok()) return
      ue=state%x(dofs)
      do k=1,5
        do j=1,5
          do i=1,5
            call hex20_geometry(coords,gp(i),gp(j),gp(k),shape,grad,det,status)
            if (.not. status%is_ok()) return
            call build_operators(shape,grad,op,testop);raw=matmul(op,ue)
            call unpack_point_state(material,raw,ps)
            call evaluate_material_point(material,ps,load,response,tangent,status)
            if (.not. status%is_ok()) return
            x=matmul(coords,shape);w=gw(i)*gw(j)*gw(k)*det
            n=n+1;weights(n)=w;emagnitude=sqrt(sum(ps%electric_field**2));efields(n)=emagnitude
            vol=vol+w;if (emagnitude>1.0e6_dp) highvol=highvol+w
            disp=0.0_dp
            do node=1,20
              disp=disp+shape(node)*ue(7*(node-1)+1:7*(node-1)+3)
            end do
            vm=sqrt(0.5_dp*((response%stress(1)-response%stress(2))**2+ &
              (response%stress(2)-response%stress(3))**2+(response%stress(3)-response%stress(1))**2)+ &
              3.0_dp*sum(response%stress(4:6)**2))
            extra(3)=max(extra(3),ps%temperature-material%reference_temperature)
            extra(4)=max(extra(4),abs(response%n/material%reference_n-1.0_dp))
            extra(5)=max(extra(5),abs(response%p/material%reference_p-1.0_dp))
            extra(6)=max(extra(6),sqrt(sum(disp**2)))
            extra(7)=min(extra(7),response%mobility_n(1,1)/material%mobility_n0(1,1))
            extra(8)=max(extra(8),vm)
            ! Full weighted 5^3 samples remain available for independent plots/comparisons.
            write(unit_id,'(*(es18.10e3,:,","))') x,w,ps%electric_field,response%joule_heat, &
              ps%temperature-material%reference_temperature,response%n/material%reference_n, &
              response%p/material%reference_p,disp,response%mobility_n(1,1)/material%mobility_n0(1,1),response%stress
          end do
        end do
      end do
    end do
    close(unit_id)
    extra(1)=weighted_percentile(efields,weights,0.95_dp);extra(2)=highvol/vol
  end subroutine extra_fields
end module f1_case_runner_m

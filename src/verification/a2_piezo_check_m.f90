module a2_piezo_check_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use field_layout_m, only : NDOF_NODE,F_U1,F_U2,F_U3,F_PHI,F_ETA_N,F_ETA_P,F_THETA,global_dof
  use model_types_m, only : MaterialType,StateType,LoadStateType,NonlinearOptionsType
  use mesh_types_m, only : MeshType
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  use hex20_shape_m, only : hex20_natural_coordinates,hex20_shape,hex20_geometry
  use mesh_checks_m, only : validate_mesh
  use sparse_pattern_m, only : SparseAssemblyContextType,initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use newton_solver_m, only : NewtonReportType,solve_single_load_newton
  use independent_flux_audit_m, only : independent_gauss_rule
  use material_point_m, only : MaterialPointStateType,MaterialPointResponseType,MaterialPointTangentType, &
    evaluate_material_point
  implicit none
  private
  real(dp), parameter :: HALF_WIDTH=10.0e-6_dp,RADIUS=3.0e-6_dp,THICKNESS=1.0e-6_dp
  real(dp), parameter :: E0=1.0e4_dp,GAMMA0=1.0e-6_dp,C44=1.0e11_dp,EPSILON=1.0e-10_dp,E15_REF=1.0_dp
  real(dp), parameter :: NREF=1.0e21_dp,MU_N=0.02_dp,MU_P=0.005_dp
  real(dp), parameter :: ZERO_W=E15_REF*E0/C44*HALF_WIDTH,ZERO_GW=E15_REF*E0/C44
  real(dp), parameter :: ZERO_PHI=E0*HALF_WIDTH,ZERO_GP=E0,ZERO_TOL=1.0e-9_dp
  integer, parameter :: NERR=15,NT(3)=[16,24,32],NR(3)=[2,3,4]
  character(len=40), parameter :: ERROR_NAMES(NERR)=[character(len=40):: &
    'w_L2','phi_L2','grad_w_L2','grad_phi_L2','stress_scaled_L2','D_scaled_L2','Q_scaled_L2', &
    'wall_traction_RMS','wall_D_RMS','wall_Jn_RMS','wall_Jp_RMS','density_n_max','density_p_max', &
    'wall_Q_peak_error','thickness_shear_traction_RMS_Pa']
  public :: verify_a2_piezo,export_a2p_matrices
contains
  subroutine export_a2p_matrices(output_dir,status)
    use model_types_m, only : DiagnosticsType
    use dirichlet_m, only : project_dirichlet_state,build_fixed_dof_mask
    use global_assembly_m, only : GlobalAssemblyRequestType,GlobalSystemType,assemble_global_system
    use scaling_m, only : build_variable_scaling,equilibrate_csr_system
    use pardiso_solver_m, only : PardisoDiagnosticsType,solve_csr_pardiso
    character(len=*), intent(in) :: output_dir
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: bc(:)
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: system
    type(DiagnosticsType) :: diagnostics
    type(PardisoDiagnosticsType) :: linear_diagnostics
    type(SolverStatusType) :: linear_status
    real(dp), allocatable :: variable_scales(:),row_scales(:),scaled_values(:),scaled_rhs(:),solution(:)
    logical, allocatable :: fixed(:)
    integer, allocatable :: fixed_int(:)
    real(dp) :: scales(7),fa,minimum_det,volume
    integer :: level,gauge,anchor,node,n,unit_id,io,gauge_dof
    character(len=48) :: label
    call clear_status(status)
    if (len_trim(output_dir)==0) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'A2PMAT requires a new existing export directory.');return
    end if
    write(*,'(a)') 'verification_id=A2PMAT_mechanical_zero_piezo_matrix_diagnostic'
    write(*,'(a)') 'scope=matrix_export_completeness_not_A2P_acceptance'
    call define_material(0.0_dp,materials(1))
    load%voltage_factor=1.0_dp;load%joule_factor=1.0_dp
    scales=[ZERO_W,ZERO_W,GAMMA0*HALF_WIDTH,ZERO_PHI,0.1_dp,0.1_dp,1.0_dp]
    request%need_jacobian=.true.
    do gauge=0,2
      do level=1,3
        write(label,'("gauge",i0,"_L",i0)') gauge,level-1
        call build_hole(NT(level),NR(level),mesh)
        call validate_mesh(mesh,minimum_det,volume,status)
        if (.not. status%is_ok()) return
        call define_bc(mesh,2,0.0_dp,bc,anchor,fa)
        if (gauge==0) then
          ! Deliberately singular diagnostic: remove only the phi gauge; never solve it.
          bc=[bc(:5),bc(7:)]
          anchor=0
        else if (gauge==1) then
          do node=1,mesh%number_of_nodes()
            if (abs(maxval(abs(mesh%coordinates(1:2,node)))-HALF_WIDTH)<1.0e-15_dp) then
              anchor=node;exit
            end if
          end do
          bc(6)%node_ids=[anchor]
        end if
        call initialize_sparse_context(mesh,context,status)
        if (.not. status%is_ok()) return
        n=context%number_of_equations
        if (allocated(state%x)) deallocate(state%x)
        allocate(state%x(n));state%x=0.0_dp
        call project_dirichlet_state(state,bc,load,mesh%number_of_nodes(),status)
        if (.not. status%is_ok()) return
        call build_fixed_dof_mask(bc,mesh%number_of_nodes(),fixed,status)
        if (.not. status%is_ok()) return
        call build_variable_scaling(mesh%number_of_nodes(),scales,variable_scales,status)
        if (.not. status%is_ok()) return
        call assemble_global_system(mesh,materials,state,load,context,bc,request,system,diagnostics,status)
        if (.not. status%is_ok()) return
        call equilibrate_csr_system(context,system%jacobian_values,-system%residual,variable_scales, &
          scaled_values,scaled_rhs,row_scales,status)
        if (.not. status%is_ok()) return
        if (allocated(solution)) deallocate(solution)
        allocate(solution(n));solution=0.0_dp
        call clear_status(linear_status)
        linear_diagnostics=PardisoDiagnosticsType()
        linear_diagnostics%perturbed_pivots=-1
        if (gauge>0) call solve_csr_pardiso(n,context%ia,context%ja,scaled_values,scaled_rhs, &
          solution,linear_diagnostics,linear_status)
        fixed_int=merge(1,0,fixed)
        gauge_dof=0
        if (anchor>0) gauge_dof=global_dof(anchor,F_PHI)
        open(newunit=unit_id,file=trim(output_dir)//'/A2PMAT_'//trim(label)//'.bin', &
          access='stream',form='unformatted',status='new',action='write',iostat=io)
        if (io/=0) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A2PMAT binary.');return
        end if
        ! Version 1: eleven int32 header values, then arrays in the documented order.
        write(unit_id) [1,n,context%number_of_nonzeros,mesh%number_of_nodes(),anchor, &
          gauge_dof,level-1,gauge,linear_diagnostics%perturbed_pivots, &
          linear_diagnostics%iterative_refinements,linear_diagnostics%error_code]
        write(unit_id) context%ia,context%ja,system%jacobian_values,scaled_values,-system%residual,scaled_rhs
        write(unit_id) variable_scales,row_scales,solution,fixed_int,mesh%coordinates
        close(unit_id)
        write(*,'(a,i0)') 'a2pmat_'//trim(label)//'_equations=',n
        write(*,'(a,i0)') 'a2pmat_'//trim(label)//'_stored_entries=',context%number_of_nonzeros
        write(*,'(a,i0)') 'a2pmat_'//trim(label)//'_gauge_dof=',gauge_dof
        write(*,'(a,i0)') 'a2pmat_'//trim(label)//'_perturbed_pivots=',linear_diagnostics%perturbed_pivots
        write(*,'(a,i0)') 'a2pmat_'//trim(label)//'_linear_error=',linear_diagnostics%error_code
        if (.not. linear_status%is_ok()) then
          status=linear_status;return
        end if
      end do
    end do
    write(*,'(a)') 'matrix_export_result=COMPLETE'
    write(*,'(a)') 'verification_result=PASS'
  end subroutine export_a2p_matrices

  subroutine verify_a2_piezo(output_dir,status)
    character(len=*), intent(in) :: output_dir
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(NonlinearOptionsType) :: options
    type(NewtonReportType) :: report
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: bc(:)
    real(dp), parameter :: signs(3)=[1.0_dp,-1.0_dp,0.0_dp]
    character(len=5), parameter :: sign_names(3)=[character(len=5)::'epos','eneg','ezero']
    character(len=10), parameter :: modes(2)=[character(len=10)::'electric','mechanical']
    real(dp) :: errors(NERR,3),scales(7),minimum_det,volume,e15,a,b,fa
    integer :: mode,sign_id,level,k,anchor,unit_summary,io,nnode
    character(len=48) :: label
    logical :: passed,case_pass,saving,gauss_pass
    call clear_status(status)
    write(*,'(a)') 'verification_id=A2P_antiplane_complementary_piezo'
    write(*,'(a)') 'model=synthetic_constrained_antiplane_nonlinear_neutral_carriers_isothermal_JE'
    write(*,'(a)') 'scope=u3_phi_solved_complementary_natural_BC_not_free_3D_plate_not_thermal_balance'
    write(*,'(a)') 'thickness_constraints=u1_u2_fixed_can_support_opposing_tangential_surface_reactions'
    write(*,'(a)') 'mechanical_contact=D_normal_zero_is_not_J_normal_zero_outer_eta_contacts_exchange_carriers'
    write(*,'(a)') 'error_quadrature=5x5x5_and_5x5_faces'
    options%max_newton=30;options%residual_tolerance=1.0e-11_dp;options%update_tolerance=1.0e-11_dp
    load%voltage_factor=1.0_dp;load%joule_factor=1.0_dp
    passed=.true.;saving=len_trim(output_dir)>0
    if (saving) then
      open(newunit=unit_summary,file=trim(output_dir)//'/A2P_summary.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A2P summary.');return
      end if
      write(unit_summary,'(a)',advance='no') 'case,level,nodes,maximum_perturbed_pivots,converged,residual,update'
      do k=1,NERR
        write(unit_summary,'(a)',advance='no') ','//trim(ERROR_NAMES(k))
      end do
      write(unit_summary,*)
    end if
    do mode=1,2
      do sign_id=1,3
        e15=signs(sign_id)*E15_REF
        call define_material(e15,materials(1))
        call exact_coefficients(mode,e15,a,b)
        scales=[ZERO_W,ZERO_W,max(abs(a)*HALF_WIDTH,ZERO_W), &
          max(abs(b)*HALF_WIDTH,ZERO_PHI),0.1_dp,0.1_dp,1.0_dp]
        case_pass=.true.
        do level=1,3
          write(label,'(a,"_",a,"_L",i0)') trim(modes(mode)),trim(sign_names(sign_id)),level-1
          call build_hole(NT(level),NR(level),mesh)
          call validate_mesh(mesh,minimum_det,volume,status)
          if (.not. status%is_ok()) return
          call define_bc(mesh,mode,e15,bc,anchor,fa)
          call initialize_sparse_context(mesh,context,status)
          if (.not. status%is_ok()) return
          nnode=mesh%number_of_nodes()
          if (allocated(state%x)) deallocate(state%x)
          allocate(state%x(NDOF_NODE*nnode));state%x=0.0_dp
          call solve_single_load_newton(mesh,materials,load,context,bc,scales,options,state,report,status)
          if (.not. status%is_ok()) return
          call measure(mesh,state,materials(1),load,mode,e15,fa,NR(level), &
            trim(output_dir),trim(label),errors(:,level),status)
          if (.not. status%is_ok()) return
          call audit_original_gauss_rows(mesh,materials,state,load,context,bc,scales, &
            options%residual_tolerance,mode,anchor,trim(output_dir),trim(label),gauss_pass,status)
          if (.not. status%is_ok()) return
          do k=1,NERR
            call metric(label,trim(ERROR_NAMES(k)),errors(k,level))
          end do
          call metric(label,'residual',report%final_residual_rms)
          call metric(label,'update',report%final_update_rms)
          write(*,'(a,i0)') 'a2p_'//trim(label)//'_nodes=',nnode
          write(*,'(a,i0)') 'a2p_'//trim(label)//'_iterations=',report%iterations
          write(*,'(a,i0)') 'a2p_'//trim(label)//'_gauge_node=',anchor
          write(*,'(a,i0)') 'a2p_'//trim(label)//'_maximum_perturbed_pivots=',report%maximum_perturbed_pivots
          call gate_metric(label,'native_converged',report%converged)
          call gate_metric(label,'zero_pivot_gate',report%maximum_perturbed_pivots==0)
          call gate_metric(label,'density_gate',maxval(errors(12:13,level))<1.0e-8_dp)
          if (a==0.0_dp) call gate_metric(label,'zero_w_gate',maxval(errors([1,3],level))<ZERO_TOL)
          if (b==0.0_dp) call gate_metric(label,'zero_phi_gate',maxval(errors([2,4],level))<ZERO_TOL)
          if (saving) write(unit_summary,'(a,",",i0,",",i0,",",i0,",",l1,*(",",es24.16e3))') &
            trim(modes(mode))//'_'//trim(sign_names(sign_id)),level-1,nnode, &
            report%maximum_perturbed_pivots,report%converged, &
            report%final_residual_rms,report%final_update_rms,errors(:,level)
          case_pass=case_pass .and. report%converged .and. report%maximum_perturbed_pivots==0 .and. gauss_pass
          case_pass=case_pass .and. maxval(errors(12:13,level))<1.0e-8_dp
          if (a==0.0_dp) case_pass=case_pass .and. maxval(errors([1,3],level))<ZERO_TOL
          if (b==0.0_dp) case_pass=case_pass .and. maxval(errors([2,4],level))<ZERO_TOL
        end do
        case_pass=case_pass .and. all(errors(1:2,3)<1.0e-3_dp) .and. all(errors(3:7,3)<0.02_dp)
        case_pass=case_pass .and. all(errors(8:11,3)<0.01_dp) .and. errors(14,3)<0.02_dp
        label=trim(modes(mode))//'_'//trim(sign_names(sign_id))
        call gate_metric(label,'final_field_gate',all(errors(1:2,3)<1.0e-3_dp))
        call gate_metric(label,'final_gradient_gate',all(errors(3:4,3)<0.02_dp))
        call gate_metric(label,'final_stress_D_Q_gate',all(errors(5:7,3)<0.02_dp))
        call gate_metric(label,'final_wall_normal_gate',all(errors(8:11,3)<0.01_dp))
        call gate_metric(label,'final_Q_peak_gate',errors(14,3)<0.02_dp)
        if (a/=0.0_dp) then
          call gate_metric(label,'w_monotonic_gate',all(errors([1,3],2)<errors([1,3],1)) &
            .and. all(errors([1,3],3)<errors([1,3],2)))
          case_pass=case_pass .and. all(errors([1,3],2)<errors([1,3],1)) &
            .and. all(errors([1,3],3)<errors([1,3],2))
        end if
        if (b/=0.0_dp) then
          call gate_metric(label,'phi_Q_J_monotonic_gate',all(errors([2,4,7,10,11],2)<errors([2,4,7,10,11],1)) &
            .and. all(errors([2,4,7,10,11],3)<errors([2,4,7,10,11],2)))
          case_pass=case_pass .and. all(errors([2,4,7,10,11],2)<errors([2,4,7,10,11],1)) &
            .and. all(errors([2,4,7,10,11],3)<errors([2,4,7,10,11],2))
        end if
        if (mode==1) then
          call gate_metric(label,'D_monotonic_gate',all(errors([6,9],2)<errors([6,9],1)) &
            .and. all(errors([6,9],3)<errors([6,9],2)))
          case_pass=case_pass .and. all(errors([6,9],2)<errors([6,9],1)) &
            .and. all(errors([6,9],3)<errors([6,9],2))
        else
          call gate_metric(label,'stress_monotonic_gate',all(errors([5,8],2)<errors([5,8],1)) &
            .and. all(errors([5,8],3)<errors([5,8],2)))
          case_pass=case_pass .and. all(errors([5,8],2)<errors([5,8],1)) &
            .and. all(errors([5,8],3)<errors([5,8],2))
        end if
        write(*,'(a,l1)') 'a2p_'//trim(modes(mode))//'_'//trim(sign_names(sign_id))//'_passed=',case_pass
        passed=passed .and. case_pass
      end do
    end do
    if (saving) close(unit_summary)
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'A2P frozen analytic or complementary-coupling criteria failed.')
    end if
  end subroutine verify_a2_piezo

  subroutine audit_original_gauss_rows(mesh,materials,state,load,context,bc,field_scales, &
                                      tolerance,mode,anchor,directory,label,passed,status)
    use model_types_m, only : DiagnosticsType
    use dirichlet_m, only : build_fixed_dof_mask
    use global_assembly_m, only : GlobalAssemblyRequestType,GlobalSystemType,assemble_global_system
    use scaling_m, only : build_variable_scaling,equilibrate_csr_system
    use independent_flux_audit_m, only : IndependentFluxReportType,integrate_independent_fluxes
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: bc(:)
    real(dp), intent(in) :: field_scales(7),tolerance
    integer, intent(in) :: mode,anchor
    character(len=*), intent(in) :: directory,label
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: raw
    type(DiagnosticsType) :: diagnostics
    type(IndependentFluxReportType) :: flux
    real(dp), allocatable :: variable_scales(:),row_scales(:),scaled_values(:),scaled_rhs(:)
    logical, allocatable :: fixed(:)
    integer, allocatable :: rows(:)
    real(dp) :: free_max,free_rms,all_max,total_scale,total_residual,identity_error,charge_error,gauge_error
    integer :: node,row,unit_id,io,nfree
    call clear_status(status)
    passed=.false.
    request%need_jacobian=.true.
    ! Reassemble the accepted state WITHOUT any Dirichlet row replacement.
    ! In particular, retain the original Gauss equation at the phi gauge.
    call assemble_global_system(mesh,materials,state,load,context,request=request, &
      system=raw,diagnostics=diagnostics,status=status)
    if (.not. status%is_ok()) return
    call build_fixed_dof_mask(bc,mesh%number_of_nodes(),fixed,status)
    if (.not. status%is_ok()) return
    call build_variable_scaling(mesh%number_of_nodes(),field_scales,variable_scales,status)
    if (.not. status%is_ok()) return
    ! Exactly the production Newton row scale s_i=max_j |J_ij|*xscale_j.
    ! These nonzero tangent scales remain meaningful at zero electric fields.
    call equilibrate_csr_system(context,raw%jacobian_values,-raw%residual,variable_scales, &
      scaled_values,scaled_rhs,row_scales,status)
    if (.not. status%is_ok()) return
    rows=[(global_dof(node,F_PHI),node=1,mesh%number_of_nodes())]
    nfree=count(.not. fixed(rows))
    if (nfree==0) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'A2P Gauss audit has no free potential rows.');return
    end if
    free_max=maxval(abs(scaled_rhs(rows)),mask=.not. fixed(rows))
    free_rms=sqrt(sum(pack(scaled_rhs(rows)**2,.not. fixed(rows)))/real(nfree,dp))
    all_max=maxval(abs(scaled_rhs(rows)))
    total_scale=sum(row_scales(rows));total_residual=sum(raw%residual(rows))
    call integrate_independent_fluxes(mesh,materials,state,load,5,flux,status)
    if (.not. status%is_ok()) return
    ! Partition of unity gives sum(original Gauss residuals)=integral(rho).
    ! Compare production assembly against independent order-five volume integration.
    identity_error=abs(total_residual-flux%charge)/total_scale
    charge_error=abs(flux%charge)/total_scale
    gauge_error=0.0_dp
    if (mode==2) gauge_error=abs(scaled_rhs(global_dof(anchor,F_PHI)))
    passed=free_max<=tolerance .and. identity_error<=tolerance
    if (mode==2) passed=passed .and. all_max<=tolerance .and. charge_error<=tolerance
    call metric(label,'original_Gauss_free_max_scaled',free_max)
    call metric(label,'original_Gauss_free_RMS_scaled',free_rms)
    call metric(label,'original_Gauss_all_max_scaled',all_max)
    call metric(label,'original_Gauss_sum_C',total_residual)
    call metric(label,'original_Gauss_total_row_scale_C',total_scale)
    call metric(label,'independent_charge_C',flux%charge)
    call metric(label,'Gauss_sum_charge_identity_scaled',identity_error)
    call metric(label,'independent_charge_scaled',charge_error)
    call metric(label,'independent_exterior_D_flux_C',flux%electric_displacement_out)
    call metric(label,'independent_exterior_D_minus_charge_C',flux%electric_displacement_out-flux%charge)
    if (mode==2) then
      call metric(label,'original_Gauss_gauge_row_scaled',gauge_error)
      call gate_metric(label,'original_Gauss_all_rows_gate',all_max<=tolerance)
      call gate_metric(label,'zero_total_charge_compatibility_gate',charge_error<=tolerance)
    else
      write(*,'(a)') 'a2p_'//trim(label)//'_outer_Gauss_rows=Dirichlet_contact_reactions_not_zero_gated'
    end if
    call gate_metric(label,'original_Gauss_free_rows_gate',free_max<=tolerance)
    call gate_metric(label,'Gauss_sum_charge_identity_gate',identity_error<=tolerance)
    call gate_metric(label,'original_Gauss_compatibility_passed',passed)
    if (len_trim(directory)>0) then
      open(newunit=unit_id,file=trim(directory)//'/A2P_'//trim(label)//'_gauss_rows.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A2P Gauss rows.');return
      end if
      write(unit_id,'(a)') 'node,x_m,y_m,z_m,potential_fixed,potential_gauge,raw_Gauss_C,row_scale_C,scaled_Gauss'
      do node=1,mesh%number_of_nodes()
        row=rows(node)
        write(unit_id,'(i0,3(",",es24.16e3),2(",",l1),3(",",es24.16e3))') &
          node,mesh%coordinates(:,node),fixed(row),mode==2.and.node==anchor, &
          raw%residual(row),row_scales(row),raw%residual(row)/row_scales(row)
      end do
      close(unit_id)
      open(newunit=unit_id,file=trim(directory)//'/A2P_'//trim(label)//'_gauss_summary.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A2P Gauss summary.');return
      end if
      write(unit_id,'(a)') 'metric,value'
      write(unit_id,'(a,es24.16e3)') 'free_max_scaled,',free_max
      write(unit_id,'(a,es24.16e3)') 'free_RMS_scaled,',free_rms
      write(unit_id,'(a,es24.16e3)') 'all_max_scaled,',all_max
      write(unit_id,'(a,es24.16e3)') 'total_row_scale_C,',total_scale
      write(unit_id,'(a,es24.16e3)') 'raw_Gauss_sum_C,',total_residual
      write(unit_id,'(a,es24.16e3)') 'independent_charge_C,',flux%charge
      write(unit_id,'(a,es24.16e3)') 'Gauss_sum_charge_identity_scaled,',identity_error
      write(unit_id,'(a,es24.16e3)') 'independent_charge_scaled,',charge_error
      write(unit_id,'(a,es24.16e3)') 'independent_exterior_D_flux_C,',flux%electric_displacement_out
      write(unit_id,'(a,es24.16e3)') 'independent_exterior_D_minus_charge_C,', &
        flux%electric_displacement_out-flux%charge
      if (mode==2) write(unit_id,'(a,es24.16e3)') 'gauge_row_scaled,',gauge_error
      write(unit_id,'(a,es24.16e3)') 'frozen_Newton_tolerance,',tolerance
      write(unit_id,'(a,l1)') 'compatibility_passed,',passed
      close(unit_id)
    end if
  end subroutine audit_original_gauss_rows

  subroutine define_material(e15,material)
    real(dp), intent(in) :: e15
    type(MaterialType), intent(out) :: material
    integer :: i
    material=MaterialType();material%name='A2P_SYNTHETIC_ANTIPLANE_JE'
    material%reference_n=NREF;material%reference_p=NREF;material%intrinsic_concentration=NREF
    material%frozen_transport=.true.;material%heat_source_current_work=.true.
    do i=1,6
      material%elastic(i,i)=C44
    end do
    do i=1,3
      material%permittivity(i,i)=EPSILON
      material%thermal_conductivity(i,i)=100.0_dp
      material%mobility_n0(i,i)=MU_N;material%mobility_p0(i,i)=MU_P
    end do
    material%piezoelectric(1,5)=e15;material%piezoelectric(2,4)=e15
  end subroutine define_material

  subroutine exact_coefficients(mode,e15,a,b)
    integer, intent(in) :: mode
    real(dp), intent(in) :: e15
    real(dp), intent(out) :: a,b
    if (mode==1) then
      a=e15*E0/C44;b=-E0
    else
      a=GAMMA0;b=e15*GAMMA0/EPSILON
    end if
  end subroutine exact_coefficients

  subroutine harmonic(x,f,g)
    real(dp), intent(in) :: x(3)
    real(dp), intent(out) :: f,g(3)
    real(dp) :: r2
    r2=sum(x(1:2)**2)
    f=x(1)*(1.0_dp+RADIUS**2/r2)
    g=[1.0_dp+RADIUS**2*(x(2)**2-x(1)**2)/r2**2, &
      -2.0_dp*RADIUS**2*x(1)*x(2)/r2**2,0.0_dp]
  end subroutine harmonic

  subroutine exact_fields(x,mode,e15,fa,w,phi,gw,gphi)
    real(dp), intent(in) :: x(3),e15,fa
    integer, intent(in) :: mode
    real(dp), intent(out) :: w,phi,gw(3),gphi(3)
    real(dp) :: a,b,f,g(3)
    call exact_coefficients(mode,e15,a,b);call harmonic(x,f,g)
    w=a*f;phi=b*f
    if (mode==1) w=a*(f-fa)
    if (mode==2) phi=b*(f-fa)
    gw=a*g;gphi=b*g
  end subroutine exact_fields

  subroutine define_bc(mesh,mode,e15,bc,anchor,fa)
    type(MeshType), intent(in) :: mesh
    integer, intent(in) :: mode
    real(dp), intent(in) :: e15
    type(DirichletBCType), allocatable, intent(out) :: bc(:)
    integer, intent(out) :: anchor
    real(dp), intent(out) :: fa
    integer, allocatable :: nodes(:),outer(:)
    integer :: i,j,driven,free_field,old_anchor
    real(dp) :: g(3),w,phi,gw(3),gphi(3)
    nodes=[(i,i=1,mesh%number_of_nodes())]
    outer=pack(nodes,abs(maxval(abs(mesh%coordinates(1:2,:)),dim=1)-HALF_WIDTH)<1.0e-15_dp)
    anchor=outer(1)
    if (mode==2) then
      ! Same physical gauge, placed symmetrically in z. No material or solver change.
      old_anchor=anchor
      do i=1,size(outer)
        j=outer(i)
        if (maxval(abs(mesh%coordinates(1:2,j)-mesh%coordinates(1:2,old_anchor)))<1.0e-15_dp &
          .and. abs(mesh%coordinates(3,j))<1.0e-15_dp) then
          anchor=j;exit
        end if
      end do
    end if
    call harmonic(mesh%coordinates(:,anchor),fa,g)
    allocate(bc(6+size(outer)))
    bc(1)%field_id=F_U1;bc(1)%node_ids=nodes
    bc(2)%field_id=F_U2;bc(2)%node_ids=nodes
    bc(3)%field_id=F_THETA;bc(3)%node_ids=nodes
    bc(4)%field_id=F_ETA_N;bc(4)%node_ids=outer
    bc(5)%field_id=F_ETA_P;bc(5)%node_ids=outer
    driven=F_PHI;free_field=F_U3
    if (mode==2) then
      driven=F_U3;free_field=F_PHI
    end if
    bc(6)%field_id=free_field;bc(6)%node_ids=[anchor]
    do j=1,size(outer)
      call exact_fields(mesh%coordinates(:,outer(j)),mode,e15,fa,w,phi,gw,gphi)
      bc(6+j)%field_id=driven;bc(6+j)%node_ids=[outer(j)]
      if (mode==1) bc(6+j)%value_end=phi
      if (mode==2) bc(6+j)%value_end=w
    end do
  end subroutine define_bc

  subroutine build_hole(ntheta,nradial,mesh)
    integer, intent(in) :: ntheta,nradial
    type(MeshType), intent(out) :: mesh
    integer, allocatable :: ids(:,:,:)
    real(dp), allocatable :: coordinates(:,:)
    real(dp) :: natural(3,20),angle,outer,rho,r,pi
    integer :: er,et,j,ir,it,iz,node,element
    pi=acos(-1.0_dp);call hex20_natural_coordinates(natural)
    allocate(ids(0:2*nradial,0:2*ntheta-1,0:2),coordinates(3,20*nradial*ntheta))
    allocate(mesh%connectivity(20,nradial*ntheta),mesh%material_id(nradial*ntheta))
    ids=0;node=0;element=0;mesh%material_id=1
    do et=0,ntheta-1
      do er=0,nradial-1
        element=element+1
        do j=1,20
          ir=2*er+nint(natural(1,j))+1
          it=modulo(2*et+nint(natural(2,j))+1,2*ntheta);iz=nint(natural(3,j))+1
          if (ids(ir,it,iz)==0) then
            node=node+1;ids(ir,it,iz)=node
            angle=pi*real(it,dp)/real(ntheta,dp)
            outer=HALF_WIDTH/max(abs(cos(angle)),abs(sin(angle)))
            rho=(real(ir,dp)/real(2*nradial,dp))**2
            if (modulo(ir,2)==1) rho=0.5_dp*((real(ir-1,dp)/real(2*nradial,dp))**2+ &
              (real(ir+1,dp)/real(2*nradial,dp))**2)
            r=RADIUS+rho*(outer-RADIUS)
            coordinates(:,node)=[r*cos(angle),r*sin(angle),THICKNESS*0.5_dp*real(iz-1,dp)]
          end if
          mesh%connectivity(j,element)=ids(ir,it,iz)
        end do
      end do
    end do
    mesh%coordinates=coordinates(:,1:node)
  end subroutine build_hole

  subroutine measure(mesh,state,material,load,mode,e15,fa,nradial,directory,label,errors,status)
    type(MeshType), intent(in) :: mesh
    type(StateType), intent(in) :: state
    type(MaterialType), intent(in) :: material
    type(LoadStateType), intent(in) :: load
    integer, intent(in) :: mode,nradial
    real(dp), intent(in) :: e15,fa
    character(len=*), intent(in) :: directory,label
    real(dp), intent(out) :: errors(NERR)
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: gp(:),weights(:)
    real(dp) :: coords(3,20),values(7,20),shape(20),gradient(3,20),dn(3,20),jac(3,3)
    real(dp) :: x(3),det,weight,volume,area,normal(3),sumwall(4),wallarea,peak_q,caparea,capstress
    real(dp) :: a,b,we,pe,gwe(3),gpe(3),wn,pn,gwn(3),gpn(3),sex(6),dex(3),qex,sigma
    real(dp) :: errsum(7),exactsum(4),sref,dref,jnref,jpref,qref,flux(4),t,angle
    integer :: element,node,i,j,k,f,unit_nodes,unit_wall,io,side
    logical :: saving
    type(MaterialPointResponseType) :: response
    call independent_gauss_rule(5,gp,weights,status)
    if (.not. status%is_ok()) return
    call exact_coefficients(mode,e15,a,b)
    sigma=ELEMENTARY_CHARGE*NREF*(MU_N+MU_P)
    sref=abs(C44*a+e15*b);dref=abs(e15*a-EPSILON*b)
    if (mode==1) sref=C44*GAMMA0+E15_REF*E0
    if (mode==2) dref=(EPSILON+E15_REF**2/C44)*E0
    jnref=ELEMENTARY_CHARGE*NREF*MU_N*abs(b);jpref=ELEMENTARY_CHARGE*NREF*MU_P*abs(b)
    qref=sigma*b*b
    if (b==0.0_dp) then
      jnref=ELEMENTARY_CHARGE*NREF*MU_N*E0;jpref=ELEMENTARY_CHARGE*NREF*MU_P*E0;qref=sigma*E0**2
    end if
    saving=len_trim(directory)>0
    if (saving) then
      open(newunit=unit_nodes,file=trim(directory)//'/A2P_'//trim(label)//'_nodes.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A2P nodes.');return
      end if
      write(unit_nodes,'(a)') 'x_m,y_m,z_m,w_m,w_exact_m,phi_V,phi_exact_V,eta_n,eta_p'
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(:,node);call exact_fields(x,mode,e15,fa,we,pe,gwe,gpe)
        write(unit_nodes,'(*(es24.16e3,:,","))') x,state%x(global_dof(node,F_U3)),we, &
          state%x(global_dof(node,F_PHI)),pe,state%x(global_dof(node,F_ETA_N)),state%x(global_dof(node,F_ETA_P))
      end do
      close(unit_nodes)
      open(newunit=unit_wall,file=trim(directory)//'/A2P_'//trim(label)//'_wall.csv', &
        status='new',action='write',iostat=io)
      if (io/=0) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A2P wall.');return
      end if
      write(unit_wall,'(a)') 'angle_rad,x_m,y_m,tz_scaled,Dn_scaled,Jn_scaled,Jp_scaled,Q_W_m3,Q_exact_W_m3'
    end if
    volume=0.0_dp;errsum=0.0_dp;exactsum=0.0_dp;sumwall=0.0_dp;wallarea=0.0_dp
    caparea=0.0_dp;capstress=0.0_dp;peak_q=0.0_dp
    do element=1,mesh%number_of_elements()
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      do node=1,20
        do f=1,7
          values(f,node)=state%x(global_dof(mesh%connectivity(node,element),f))
        end do
      end do
      do k=1,5
        do j=1,5
          do i=1,5
            call hex20_geometry(coords,gp(i),gp(j),gp(k),shape,gradient,det,status)
            if (.not. status%is_ok()) return
            weight=weights(i)*weights(j)*weights(k)*det;volume=volume+weight
            call point_values()
            if (.not. status%is_ok()) return
            errsum=errsum+weight*[(wn-we)**2,(pn-pe)**2,sum((gwn-gwe)**2),sum((gpn-gpe)**2), &
              sum((response%stress-sex)**2)/sref**2,sum((response%electric_displacement-dex)**2)/dref**2, &
              ((response%joule_heat-qex)/qref)**2]
            exactsum=exactsum+weight*[we**2,pe**2,sum(gwe**2),sum(gpe**2)]
          end do
        end do
      end do
      ! Direct surface quadrature: constraints on u1/u2 can carry these shear tractions.
      do side=-1,1,2
        do j=1,5
          do i=1,5
            call hex20_geometry(coords,gp(i),gp(j),real(side,dp),shape,gradient,det,status)
            if (.not. status%is_ok()) return
            call point_values()
            if (.not. status%is_ok()) return
            call hex20_shape(gp(i),gp(j),real(side,dp),shape,dn)
            jac=matmul(coords,transpose(dn));area=sqrt(sum(cross(jac(:,1),jac(:,2))**2))
            weight=weights(i)*weights(j)*area
            caparea=caparea+weight;capstress=capstress+weight*sum(response%stress(4:5)**2)
          end do
        end do
      end do
      if (modulo(element-1,nradial)/=0) cycle
      do k=1,5
        do j=1,5
          call wall_point(gp(j),gp(k))
          if (.not. status%is_ok()) return
          weight=weights(j)*weights(k)*area
          sumwall=sumwall+weight*flux**2;wallarea=wallarea+weight
        end do
      end do
      do j=0,20
        t=-1.0_dp+0.1_dp*real(j,dp);call wall_point(t,0.0_dp)
        if (.not. status%is_ok()) return
        peak_q=max(peak_q,response%joule_heat/qref);angle=atan2(x(2),x(1))
        if (saving) write(unit_wall,'(*(es24.16e3,:,","))') angle,x(1:2),flux,response%joule_heat,qex
      end do
    end do
    if (saving) close(unit_wall)
    if (a==0.0_dp) then
      exactsum(1)=volume*ZERO_W**2;exactsum(3)=volume*ZERO_GW**2
    end if
    if (b==0.0_dp) then
      exactsum(2)=volume*ZERO_PHI**2;exactsum(4)=volume*ZERO_GP**2
    end if
    errors(1:4)=sqrt(errsum(1:4)/exactsum)
    errors(5:7)=sqrt(errsum(5:7)/volume);errors(8:11)=sqrt(sumwall/wallarea)
    errors(12)=maxval(abs(exp(state%x(F_ETA_N::NDOF_NODE))-1.0_dp))
    errors(13)=maxval(abs(exp(state%x(F_ETA_P::NDOF_NODE))-1.0_dp))
    errors(14)=abs(peak_q/4.0_dp-1.0_dp)
    if (b==0.0_dp) errors(14)=abs(peak_q)
    errors(15)=sqrt(capstress/caparea)
  contains
    subroutine point_values()
      type(MaterialPointStateType) :: point
      type(MaterialPointTangentType) :: tangent
      real(dp) :: gu(3,3)
      integer :: ii
      x=matmul(coords,shape);call exact_fields(x,mode,e15,fa,we,pe,gwe,gpe)
      wn=dot_product(shape,values(F_U3,:));pn=dot_product(shape,values(F_PHI,:))
      gwn=matmul(gradient,values(F_U3,:));gpn=matmul(gradient,values(F_PHI,:))
      do ii=1,3
        gu(:,ii)=matmul(gradient,values(ii,:))
      end do
      point=MaterialPointStateType()
      point%strain=[gu(1,1),gu(2,2),gu(3,3),gu(3,2)+gu(2,3),gu(3,1)+gu(1,3),gu(2,1)+gu(1,2)]
      point%electric_field=-gpn;point%temperature=material%reference_temperature+dot_product(shape,values(F_THETA,:))
      point%temperature_gradient=matmul(gradient,values(F_THETA,:))
      point%eta_n=dot_product(shape,values(F_ETA_N,:));point%eta_p=dot_product(shape,values(F_ETA_P,:))
      point%eta_n_gradient=matmul(gradient,values(F_ETA_N,:));point%eta_p_gradient=matmul(gradient,values(F_ETA_P,:))
      call evaluate_material_point(material,point,load,response,tangent,status)
      sex=0.0_dp;sex(4)=C44*gwe(2)+e15*gpe(2);sex(5)=C44*gwe(1)+e15*gpe(1)
      dex=e15*gwe-EPSILON*gpe;qex=sigma*sum(gpe**2)
    end subroutine point_values

    subroutine wall_point(t,z)
      real(dp), intent(in) :: t,z
      call hex20_geometry(coords,-1.0_dp,t,z,shape,gradient,det,status)
      if (.not. status%is_ok()) return
      call point_values()
      if (.not. status%is_ok()) return
      call hex20_shape(-1.0_dp,t,z,shape,dn)
      jac=matmul(coords,transpose(dn));normal=cross(jac(:,2),jac(:,3))
      area=sqrt(sum(normal**2));normal=-normal/area
      flux=[dot_product([response%stress(5),response%stress(4),response%stress(3)],normal)/sref, &
        dot_product(response%electric_displacement,normal)/dref, &
        dot_product(response%current_n,normal)/jnref,dot_product(response%current_p,normal)/jpref]
    end subroutine wall_point
  end subroutine measure

  pure function cross(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function cross

  subroutine metric(label,name,value)
    character(len=*), intent(in) :: label,name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') 'a2p_'//trim(label)//'_'//trim(name)//'=',value
  end subroutine metric

  subroutine gate_metric(label,name,passed)
    character(len=*), intent(in) :: label,name
    logical, intent(in) :: passed
    write(*,'(a,l1)') 'a2p_'//trim(label)//'_'//trim(name)//'=',passed
  end subroutine gate_metric
end module a2_piezo_check_m

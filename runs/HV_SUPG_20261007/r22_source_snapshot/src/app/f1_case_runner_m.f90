module f1_case_runner_m
  use, intrinsic :: iso_fortran_env, only : output_unit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use field_layout_m, only : NDOF_NODE,F_U1,F_U2,F_U3,F_PHI,F_ETA_N,F_ETA_P, &
    F_THETA,global_dof
  use input_schema_m, only : CaseConfigType
  use model_types_m, only : MaterialType,StateType,LoadStateType,DiagnosticsType
  use mesh_types_m, only : MeshType
  use native_input_m, only : read_native_case
  use mesh_checks_m, only : validate_mesh
  use zno_candidate_m, only : define_zno_candidate
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_INPUT_ERROR,STATUS_VERIFICATION_FAILED, &
    STATUS_NONLINEAR_SOLVER_FAILED
  use dirichlet_m, only : DirichletBCType,build_fixed_dof_mask
  use sparse_pattern_m, only : SparseAssemblyContextType,initialize_sparse_context
  use newton_solver_m, only : NewtonReportType,solve_single_load_newton
  use temperature_control_m, only : solve_temperature_controlled_newton
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
  subroutine run_f1_case(case_path,directory,status,inverse_step,inverse_count)
    character(len=*), intent(in) :: case_path,directory
    type(SolverStatusType), intent(out) :: status
    real(dp), optional, intent(in) :: inverse_step
    integer, optional, intent(in) :: inverse_count
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
    real(dp) :: heat_error,jn_error,jp_error,quad_error,power_scale,current_scale,extra(8),domain_fraction
    real(dp) :: external_heat_error,external_jn_error,external_jp_error,symmetry_heat_error,symmetry_jn_error,symmetry_jp_error
    integer :: unit_id,io,stage,point,total_iterations,inverse_node
    logical :: accepted
    call clear_status(status)
    domain_fraction=1.0_dp
    inverse_node=0
    if(present(inverse_step).neqv.present(inverse_count)) then
      call set_status(status,STATUS_INPUT_ERROR,'Inverse diagnostics require both temperature step and point count.');return
    end if
    if(present(inverse_step)) then
      if(.not.ieee_is_finite(inverse_step).or.inverse_step<=0.or.inverse_count<1.or.inverse_count>100) then
        call set_status(status,STATUS_INPUT_ERROR,'Invalid inverse diagnostic step or count.');return
      end if
    end if
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
    if(config%half_y_symmetry) then
      domain_fraction=0.5_dp
      write(unit_id,'(a)') 'domain_symmetry=Y_HALF'
      write(unit_id,'(a)') 'integral_scope=ACTUAL_HALF_DOMAIN_NO_AUTOMATIC_DOUBLING'
    end if
    flush(unit_id)
    call validate_mesh(mesh,minimum_det,volume,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call define_zno_candidate(model,materials(1))
    materials(1)%transport_supg=config%transport_supg
    materials(1)%transport_log_weighted=config%transport_log_weighted
    materials(1)%supg_scale=config%supg_scale
    materials(1)%supg_directional_tau=config%supg_directional_tau
    if(config%transport_log_weighted) then
      write(unit_id,'(a)') 'transport_discretization=HEX20_LOG_WEIGHTED_EXPERIMENTAL'
      write(unit_id,'(a)') 'supg_parameters_applicable=F'
      write(unit_id,'(a,es24.16e3)') 'supg_scale=',config%supg_scale
      write(unit_id,'(a,l1)') 'supg_directional_tau=',config%supg_directional_tau
      write(unit_id,'(a)') 'algebraic_carrier_balance=WEIGHTED_FREE_ROWS_NOT_PHYSICAL_CURRENT'
    else if(config%transport_supg) then
      write(unit_id,'(a)') 'transport_discretization=HEX20_CONSISTENT_SUPG_EXPERIMENTAL'
      write(unit_id,'(a,es24.16e3)') 'supg_scale=',config%supg_scale
      write(unit_id,'(a,l1)') 'supg_directional_tau=',config%supg_directional_tau
    end if
    config%continuation_options%heat_balance_reference_power=ELEMENTARY_CHARGE*( &
      materials(1)%reference_n*materials(1)%mobility_n0(1,1)+ &
      materials(1)%reference_p*materials(1)%mobility_p0(1,1))*config%target_voltage**2*1.0e-6_dp
    if(config%half_y_symmetry) config%continuation_options%heat_balance_reference_power= &
      domain_fraction*config%continuation_options%heat_balance_reference_power
    write(unit_id,'(a,es24.16e3)') 'zero_source_heat_reference_W=', &
      config%continuation_options%heat_balance_reference_power
    if(config%transport_supg.or.config%transport_log_weighted) then
      write(output_unit,'(a)') 'f1_setup=BOUNDARIES';flush(output_unit)
    end if
    call define_boundaries(mesh,config%target_voltage,bc,status,config%half_y_symmetry)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    if(config%transport_supg.or.config%transport_log_weighted) then
      write(output_unit,'(a)') 'f1_setup=SPARSE_PATTERN';flush(output_unit)
    end if
    call initialize_sparse_context(mesh,context,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    write(output_unit,'(a,i0,a,i0)') 'f1_nodes=',mesh%number_of_nodes(),', elements=',mesh%number_of_elements()
    flush(output_unit)
    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()));state%x=0.0_dp
    if (len_trim(config%restart_nodes_file) > 0) then
      call read_restart_nodes(config%restart_nodes_file,mesh,state,status)
      if (.not. status%is_ok()) then
        close(unit_id);return
      end if
      write(unit_id,'(a)') 'restart_nodes_file='//trim(config%restart_nodes_file)
    end if
    config%nonlinear_options%verbose=.true.
    max_balance=0.0_dp;max_residual=0.0_dp;max_update=0.0_dp;total_iterations=0
    if(present(inverse_step)) then
      if(model=='LQ'.or..not.config%transport_supg.or.len_trim(config%restart_nodes_file)==0.or. &
         config%restart_joule_factor/=1.or.config%restart_voltage_factor<=0.or.config%restart_voltage_factor>=1) then
        call set_status(status,STATUS_INPUT_ERROR,'Inverse diagnostics require a nonzero, subtarget full-heat SUPG restart.')
      else
        call trace_temperature_branch()
      end if
      close(unit_id);return
    end if
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
      load%voltage_factor=config%restart_voltage_factor
      load%joule_factor=config%restart_joule_factor
      target%voltage_factor=config%target_voltage_factor
      target%joule_factor=config%target_joule_factor
      if (len_trim(config%restart_nodes_file) == 0 .and. &
          config%continuation_options%full_joule_during_voltage) then
        load%joule_factor=target%joule_factor
      end if
      if (config%continuation_options%full_joule_during_voltage) then
        write(unit_id,'(a)') 'continuation_path=FULL_JOULE_DURING_VOLTAGE'
      else
        write(unit_id,'(a)') 'continuation_path=VOLTAGE_THEN_JOULE'
      end if
      if(config%transport_supg.or.config%transport_log_weighted) then
        call solve_load_continuation(mesh,materials,context,bc,config%field_scales,config%nonlinear_options, &
          config%continuation_options,target,state,load,continuation,status,save_accepted)
      else
        call solve_load_continuation(mesh,materials,context,bc,config%field_scales,config%nonlinear_options, &
          config%continuation_options,target,state,load,continuation,status)
      end if
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
    call integrate_independent_fluxes(mesh,materials,state,load,5,flux,status,half_y=config%half_y_symmetry)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call integrate_independent_fluxes(mesh,materials,state,load,4,flux4,status,half_y=config%half_y_symmetry)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    if(config%half_y_symmetry) then
      call evaluate_hole_field_metrics(mesh,materials,state,load,config%hole_radius,0.5e-6_dp,ring,status, &
        hole_center=[0.0_dp,0.0_dp])
    else
      call evaluate_hole_field_metrics(mesh,materials,state,load,config%hole_radius,0.5e-6_dp,ring,status)
    end if
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    call extra_fields(directory,mesh,materials(1),state,load,extra,status)
    if (.not. status%is_ok()) then
      close(unit_id);return
    end if
    power_scale=max(abs(flux%joule),abs(sum(flux%heat_out)),1.7736095e-3_dp*1.0e-12_dp)
    if(config%half_y_symmetry) power_scale=max(abs(flux%joule),abs(sum(flux%heat_out)), &
      domain_fraction*1.7736095e-3_dp*1.0e-12_dp)
    current_scale=max(sum(abs(flux%electron_out(1:2))),tiny(1.0_dp))
    heat_error=abs(sum(flux%heat_out)-flux%joule)/power_scale
    jn_error=abs(sum(flux%electron_out))/current_scale
    jp_error=abs(sum(flux%hole_out))/max(sum(abs(flux%hole_out(1:2))),tiny(1.0_dp))
    quad_error=abs(flux%joule-flux4%joule)/max(abs(flux%joule),1.7736095e-15_dp)
    if(config%half_y_symmetry) quad_error=abs(flux%joule-flux4%joule)/ &
      max(abs(flux%joule),domain_fraction*1.7736095e-15_dp)
    call metric('voltage_V',abs(config%target_voltage))
    call metric('hole_radius_m',config%hole_radius)
    call metric('ring_width_m',0.5e-6_dp)
    call metric('minimum_jacobian',minimum_det)
    call metric('total_volume_m3',flux%volume)
    call metric('maximum_stage_balance',max_balance)
    call metric('maximum_accepted_residual',max_residual)
    call metric('maximum_accepted_update',max_update)
    if(config%transport_supg.or.config%transport_log_weighted) then
      call metric('final_voltage_factor',load%voltage_factor)
      call metric('final_joule_factor',load%joule_factor)
      write(unit_id,'(a)') 'spatial_acceptance=REQUIRES_INDEPENDENT_MESH_AND_STABILIZATION_STUDIES'
    end if
    call metric('joule_W',flux%joule)
    call metric('terminal_power_W',flux%terminal_power)
    call metric('heat_out_W',sum(flux%heat_out))
    call metric('heat_balance_relative',heat_error)
    call metric('electron_balance_relative',jn_error)
    call metric('hole_balance_diagnostic',jp_error)
    if(config%half_y_symmetry) then
      ! Retain the actual-half balance above, and independently require the
      ! original physical exterior and artificial symmetry cut to pass. Their
      ! errors must not cancel to grant a full-domain conservation claim.
      power_scale=max(abs(flux%joule),abs(sum(flux%heat_out)-flux%symmetry_heat_out), &
        domain_fraction*1.7736095e-15_dp)
      external_heat_error=abs(sum(flux%heat_out)-flux%symmetry_heat_out-flux%joule)/power_scale
      external_jn_error=abs(sum(flux%electron_out)-flux%symmetry_electron_out)/current_scale
      external_jp_error=abs(sum(flux%hole_out)-flux%symmetry_hole_out)/ &
        max(sum(abs(flux%hole_out(1:2))),tiny(1.0_dp))
      symmetry_heat_error=abs(flux%symmetry_heat_out)/power_scale
      symmetry_jn_error=abs(flux%symmetry_electron_out)/current_scale
      symmetry_jp_error=abs(flux%symmetry_hole_out)/max(sum(abs(flux%hole_out(1:2))),tiny(1.0_dp))
      call metric('physical_exterior_heat_out_W',sum(flux%heat_out)-flux%symmetry_heat_out)
      call metric('physical_exterior_terminal_power_W',flux%terminal_power-flux%symmetry_terminal_power)
      call metric('physical_exterior_heat_balance_relative',external_heat_error)
      call metric('physical_exterior_electron_balance_relative',external_jn_error)
      call metric('physical_exterior_hole_balance_diagnostic',external_jp_error)
      call metric('symmetry_heat_out_W',flux%symmetry_heat_out)
      call metric('symmetry_electron_out_A',flux%symmetry_electron_out)
      call metric('symmetry_hole_out_A',flux%symmetry_hole_out)
      call metric('symmetry_heat_flux_relative',symmetry_heat_error)
      call metric('symmetry_electron_flux_relative',symmetry_jn_error)
      call metric('symmetry_hole_flux_diagnostic',symmetry_jp_error)
    end if
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
    if(config%transport_supg.or.config%transport_log_weighted) accepted=accepted .and. jp_error<0.01_dp .and. &
      load%voltage_factor==1.0_dp .and. load%joule_factor==1.0_dp
    if(config%half_y_symmetry) accepted=accepted .and. &
      external_heat_error<0.01_dp.and.external_jn_error<0.01_dp.and.external_jp_error<0.01_dp.and. &
      symmetry_heat_error<0.01_dp.and.symmetry_jn_error<0.01_dp.and.symmetry_jp_error<0.01_dp
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
    subroutine trace_temperature_branch()
      type(StateType) :: base_state,previous_state
      type(LoadStateType) :: base_load,previous_load
      type(SolverStatusType) :: attempt_status
      logical, allocatable :: fixed(:)
      real(dp) :: base_theta,previous_theta,step,theta_goal,ratio,max_factor
      integer :: node,equation,trace_unit,accepted_count,attempts,descending
      logical :: have_previous,saw_increase,stopped_by_criterion
      write(unit_id,'(a)') 'solution_mode=INVERSE_TEMPERATURE_DIAGNOSTIC_ONLY'
      write(unit_id,'(a)') 'acceptance=REQUIRES_FIXED_VOLTAGE_REPLAY_FLUX_AND_MESH_AUDITS'
      load%voltage_factor=config%restart_voltage_factor;load%joule_factor=1
      call solve_single_load_newton(mesh,materials,load,context,bc,config%field_scales, &
        config%nonlinear_options,state,newton,status)
      if(.not.status%is_ok()) return
      call build_fixed_dof_mask(bc,mesh%number_of_nodes(),fixed,status)
      if(.not.status%is_ok()) return
      base_theta=-huge(1.0_dp)
      do node=1,mesh%number_of_nodes()
        equation=global_dof(node,F_THETA)
        if(fixed(equation)) cycle
        if(state%x(equation)>base_theta) then
          base_theta=state%x(equation);inverse_node=node
        end if
      end do
      if(base_theta<=0.or.inverse_node==0) then
        call set_status(status,STATUS_INPUT_ERROR,'Inverse control needs a warm, free temperature node.');return
      end if
      equation=global_dof(inverse_node,F_THETA)
      write(unit_id,'(a,i0)') 'inverse_control_node=',inverse_node
      open(newunit=trace_unit,file=trim(directory)//'/inverse_branch.csv',status='new',action='write',iostat=io)
      if(io/=0) then
        call set_status(status,STATUS_INPUT_ERROR,'Cannot create inverse branch trace.');return
      end if
      write(trace_unit,'(a)') 'point,voltage_V,theta_control_K,residual,update,newton_iterations'
      call save_accepted(state,load,1,newton%final_residual_rms,newton%final_update_rms,status)
      if(.not.status%is_ok()) then
        close(trace_unit);return
      end if
      write(trace_unit,'(i0,4(",",es24.16e3),",",i0)') 1,abs(config%target_voltage*load%voltage_factor), &
        base_theta,newton%final_residual_rms,newton%final_update_rms,newton%iterations
      flush(trace_unit)
      base_state=state;base_load=load;step=inverse_step
      accepted_count=0;attempts=0;descending=0;have_previous=.false.;saw_increase=.false.
      stopped_by_criterion=.false.
      max_factor=load%voltage_factor
      do while(accepted_count<inverse_count.and.attempts<20*inverse_count)
        attempts=attempts+1;theta_goal=base_theta+step;state=base_state;load=base_load
        if(have_previous) then
          ratio=step/(base_theta-previous_theta)
          state%x=base_state%x+ratio*(base_state%x-previous_state%x)
          load%voltage_factor=base_load%voltage_factor+ratio*(base_load%voltage_factor-previous_load%voltage_factor)
        else
          state%x(F_THETA::NDOF_NODE)=state%x(F_THETA::NDOF_NODE)*(theta_goal/base_theta)
          load%voltage_factor=base_load%voltage_factor*sqrt(theta_goal/base_theta)
          ! Lift the whole potential guess when changing the electrode factor;
          ! changing contacts alone creates a spurious field in the first layer.
          do node=1,mesh%number_of_nodes()
            ratio=(mesh%coordinates(1,node)+10e-6_dp)/(20e-6_dp)
            state%x(global_dof(node,F_PHI))=state%x(global_dof(node,F_PHI))+ &
              (load%voltage_factor-base_load%voltage_factor)*config%target_voltage*ratio
          end do
        end if
        call solve_temperature_controlled_newton(mesh,materials,load,context,bc,config%field_scales, &
          config%nonlinear_options,inverse_node,theta_goal,state,newton,attempt_status)
        if(attempt_status%is_ok()) then
          call audit_global_balances(mesh,materials,state,load,context,bc,balance,bc_error,diagnostics,attempt_status)
          if(attempt_status%is_ok()) then
            if(maxval(balance)>=1e-8_dp.or.bc_error>=1e-10_dp) &
              call set_status(attempt_status,STATUS_VERIFICATION_FAILED,'Inverse point fails original algebraic balances.')
          end if
        end if
        if(.not.attempt_status%is_ok()) then
          state=base_state;load=base_load;step=.5_dp*step
          write(output_unit,'(a,es16.8,a,a)') 'inverse_cutback_K=',step,',reason=',trim(attempt_status%message)
          flush(output_unit)
          if(step<inverse_step/64) then
            status=attempt_status;exit
          end if
          cycle
        end if
        accepted_count=accepted_count+1
        call save_accepted(state,load,accepted_count+1,newton%final_residual_rms,newton%final_update_rms,status)
        if(.not.status%is_ok()) exit
        write(trace_unit,'(i0,4(",",es24.16e3),",",i0)') accepted_count+1, &
          abs(config%target_voltage*load%voltage_factor),theta_goal,newton%final_residual_rms, &
          newton%final_update_rms,newton%iterations
        flush(trace_unit)
        if(load%voltage_factor>base_load%voltage_factor+1e-5_dp) saw_increase=.true.
        if(load%voltage_factor<base_load%voltage_factor-1e-5_dp) then
          descending=descending+1
        else
          descending=0
        end if
        max_factor=max(max_factor,load%voltage_factor)
        previous_state=base_state;previous_load=base_load;previous_theta=base_theta
        base_state=state;base_load=load;base_theta=theta_goal;have_previous=.true.
        step=min(inverse_step,step*1.5_dp)
        if(load%voltage_factor>=1) then
          stopped_by_criterion=.true.
          write(unit_id,'(a)') 'inverse_stop=TARGET_CROSSED_REQUIRES_FIXED_TARGET_REPLAY';exit
        end if
        if(descending>=3) then
          stopped_by_criterion=.true.
          if(saw_increase) then
            write(unit_id,'(a)') 'inverse_stop=SAMPLED_TURNING_BRANCH_NOT_A_DYNAMICAL_STABILITY_PROOF'
          else
            write(unit_id,'(a)') 'inverse_stop=DESCENDING_BRANCH_REQUIRES_LOWER_TEMPERATURE_BRACKET'
          end if
          exit
        end if
      end do
      if(status%is_ok().and..not.stopped_by_criterion.and.accepted_count<inverse_count) then
        write(unit_id,'(a)') 'inverse_stop=ATTEMPT_LIMIT'
        call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED,'Inverse trace reached its attempt limit; partial points retained.')
      end if
      call metric('inverse_maximum_sampled_voltage_V',abs(config%target_voltage)*max_factor)
      write(unit_id,'(a,i0)') 'inverse_accepted_increments=',accepted_count
      write(unit_id,'(a)') 'inverse_scope=ALGEBRAIC_BRANCH_ONLY_NO_COMPLETE_0_8V_CLAIM'
      close(trace_unit)
    end subroutine trace_temperature_branch

    subroutine save_accepted(saved_state,saved_load,index,residual,update,saved_status)
      type(StateType), intent(in) :: saved_state
      type(LoadStateType), intent(in) :: saved_load
      integer, intent(in) :: index
      real(dp), intent(in) :: residual,update
      type(SolverStatusType), intent(out) :: saved_status
      integer :: checkpoint_unit,checkpoint_io
      character(len=40) :: name
      write(name,'(a,i5.5)') 'accepted_',index
      call write_nodes(directory,mesh,materials(1),saved_state,saved_status,trim(name)//'_nodes.csv')
      if(.not.saved_status%is_ok()) return
      open(newunit=checkpoint_unit,file=trim(directory)//'/'//trim(name)//'_load.txt', &
        status='new',action='write',iostat=checkpoint_io)
      if(checkpoint_io/=0) then
        call set_status(saved_status,STATUS_INPUT_ERROR,'Cannot create accepted-load checkpoint.');return
      end if
      write(checkpoint_unit,'(a)') 'status=ALGEBRAIC_CHECKPOINT_NOT_FLUX_OR_MESH_ACCEPTANCE'
      write(checkpoint_unit,'(a,es24.16e3)') 'applied_voltage_V=',config%target_voltage*saved_load%voltage_factor
      write(checkpoint_unit,'(a,es24.16e3)') 'voltage_factor=',saved_load%voltage_factor
      write(checkpoint_unit,'(a,es24.16e3)') 'joule_factor=',saved_load%joule_factor
      write(checkpoint_unit,'(a,es24.16e3)') 'residual=',residual
      write(checkpoint_unit,'(a,es24.16e3)') 'update=',update
      if(inverse_node>0) then
        write(checkpoint_unit,'(a,i0)') 'temperature_control_node=',inverse_node
        write(checkpoint_unit,'(a,es24.16e3)') 'theta_control_K=',saved_state%x(global_dof(inverse_node,F_THETA))
        write(checkpoint_unit,'(a)') 'replay=REQUIRES_ORDINARY_FIXED_VOLTAGE_SOLVE'
      end if
      close(checkpoint_unit)
    end subroutine save_accepted
    subroutine metric(name,value)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      write(unit_id,'(a,es24.16e3)') trim(name)//'=',value
      write(output_unit,'(a,es24.16e3)') trim(name)//'=',value
    end subroutine metric
  end subroutine run_f1_case

  subroutine read_restart_nodes(path,mesh,state,status)
    character(len=*), intent(in) :: path
    type(MeshType), intent(in) :: mesh
    type(StateType), intent(inout) :: state
    type(SolverStatusType), intent(out) :: status
    character(len=512) :: header,io_message
    real(dp) :: x,y,z,u1,u2,u3,phi,eta_n,eta_p,theta,n_value,p_value
    real(dp) :: coordinate_scale
    integer :: unit_id,io,node,node_id

    call clear_status(status)
    open(newunit=unit_id,file=trim(path),status='old',action='read',iostat=io,iomsg=io_message)
    if (io /= 0) then
      call set_status(status,STATUS_INPUT_ERROR,'Cannot open restart nodes file: '//trim(io_message))
      return
    end if
    read(unit_id,'(a)',iostat=io,iomsg=io_message) header
    if (io /= 0 .or. index(trim(header),'node,x_m,y_m,z_m,ux_m,uy_m,uz_m,phi_V,') /= 1) then
      close(unit_id)
      call set_status(status,STATUS_INPUT_ERROR,'Restart nodes file has an incompatible header.')
      return
    end if
    do node=1,mesh%number_of_nodes()
      read(unit_id,*,iostat=io,iomsg=io_message) node_id,x,y,z,u1,u2,u3,phi, &
        eta_n,eta_p,theta,n_value,p_value
      if (io /= 0 .or. node_id /= node) then
        close(unit_id)
        call set_status(status,STATUS_INPUT_ERROR, &
          'Restart nodes file is incomplete or node identifiers are not ordered.')
        return
      end if
      coordinate_scale=max(1.0e-6_dp,maxval(abs(mesh%coordinates(:,node))))
      if (maxval(abs([x,y,z]-mesh%coordinates(:,node))) > 1.0e-12_dp*coordinate_scale) then
        close(unit_id)
        call set_status(status,STATUS_INPUT_ERROR,'Restart nodes do not match the current mesh.')
        return
      end if
      state%x(global_dof(node,F_U1))=u1
      state%x(global_dof(node,F_U2))=u2
      state%x(global_dof(node,F_U3))=u3
      state%x(global_dof(node,F_PHI))=phi
      state%x(global_dof(node,F_ETA_N))=eta_n
      state%x(global_dof(node,F_ETA_P))=eta_p
      state%x(global_dof(node,F_THETA))=theta
    end do
    read(unit_id,'(a)',iostat=io) header
    close(unit_id)
    if (io == 0) then
      call set_status(status,STATUS_INPUT_ERROR,'Restart nodes file contains extra data rows.')
    end if
  end subroutine read_restart_nodes

  subroutine define_boundaries(mesh,voltage,bc,status,half_y)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: voltage
    type(DirichletBCType), allocatable, intent(out) :: bc(:)
    type(SolverStatusType), intent(out) :: status
    logical, intent(in) :: half_y
    integer, allocatable :: left(:),right(:)
    integer :: a(1),b(1),c(1),anchors(3),nleft,nright,il,ir
    real(dp) :: xyz(3)
    real(dp), parameter :: tol=1.0e-15_dp
    integer :: i
    call clear_status(status)
    if(half_y) then
      call define_half_y_boundaries(mesh,voltage,bc,status)
      return
    end if
    ! Scalar scans avoid several 3*node_count automatic SPREAD temporaries.
    ! The geometric predicates and ascending contact-node order are unchanged.
    nleft=0;nright=0;anchors=0;a=0;b=0;c=0
    do i=1,mesh%number_of_nodes()
      xyz=mesh%coordinates(:,i)
      if(abs(xyz(1)+10e-6_dp)<tol) nleft=nleft+1
      if(abs(xyz(1)-10e-6_dp)<tol) nright=nright+1
      if(maxval(abs(xyz-[-10e-6_dp,-10e-6_dp,-0.5e-6_dp]))<tol) then
        anchors(1)=anchors(1)+1;a=i
      end if
      if(maxval(abs(xyz-[-10e-6_dp,10e-6_dp,-0.5e-6_dp]))<tol) then
        anchors(2)=anchors(2)+1;b=i
      end if
      if(maxval(abs(xyz-[-10e-6_dp,-10e-6_dp,0.5e-6_dp]))<tol) then
        anchors(3)=anchors(3)+1;c=i
      end if
    end do
    if (any(anchors/=1) .or. nleft<8 .or. nright<8) then
      call set_status(status,STATUS_INPUT_ERROR,'F1 mesh lacks prescribed contacts or rigid-body anchors.');return
    end if
    allocate(left(nleft),right(nright));il=0;ir=0
    do i=1,mesh%number_of_nodes()
      if(abs(mesh%coordinates(1,i)+10e-6_dp)<tol) then
        il=il+1;left(il)=i
      end if
      if(abs(mesh%coordinates(1,i)-10e-6_dp)<tol) then
        ir=ir+1;right(ir)=i
      end if
    end do
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

  subroutine define_half_y_boundaries(mesh,voltage,bc,status)
    ! Reflection plane y=0: uy=0, natural tangential tractions and scalar fluxes.
    ! Only the remaining x/z translations and y rotation are removed by A/B.
    ! Retaining the full-domain anchors would overconstrain the half problem.
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: voltage
    type(DirichletBCType), allocatable, intent(out) :: bc(:)
    type(SolverStatusType), intent(out) :: status
    integer, allocatable :: left(:),right(:),symmetry(:)
    integer :: i,a,b,na,nb,nleft,nright,nsym,il,ir,isym
    real(dp) :: xyz(3),lo(3),hi(3)
    real(dp), parameter :: tol=1.0e-15_dp
    call clear_status(status)
    lo=minval(mesh%coordinates,dim=2);hi=maxval(mesh%coordinates,dim=2)
    if(maxval(abs(lo-[-10e-6_dp,0.0_dp,-0.5e-6_dp]))>tol.or. &
       maxval(abs(hi-[10e-6_dp,10e-6_dp,0.5e-6_dp]))>tol) then
      call set_status(status,STATUS_INPUT_ERROR,'Y-half mesh must occupy x=[-10,10], y=[0,10], z=[-0.5,0.5] um.')
      return
    end if
    nleft=0;nright=0;nsym=0;na=0;nb=0;a=0;b=0
    do i=1,mesh%number_of_nodes()
      xyz=mesh%coordinates(:,i)
      if(abs(xyz(1)+10e-6_dp)<tol) nleft=nleft+1
      if(abs(xyz(1)-10e-6_dp)<tol) nright=nright+1
      if(abs(xyz(2))<tol) nsym=nsym+1
      if(maxval(abs(xyz-[-10e-6_dp,0.0_dp,-0.5e-6_dp]))<tol) then
        a=i;na=na+1
      end if
      if(maxval(abs(xyz-[10e-6_dp,0.0_dp,-0.5e-6_dp]))<tol) then
        b=i;nb=nb+1
      end if
    end do
    if(na/=1.or.nb/=1.or.nleft<8.or.nright<8.or.nsym<8) then
      call set_status(status,STATUS_INPUT_ERROR,'Y-half mesh lacks contacts, symmetry nodes or in-plane rigid anchors.')
      return
    end if
    allocate(left(nleft),right(nright),symmetry(nsym));il=0;ir=0;isym=0
    do i=1,mesh%number_of_nodes()
      xyz=mesh%coordinates(:,i)
      if(abs(xyz(1)+10e-6_dp)<tol) then
        il=il+1;left(il)=i
      end if
      if(abs(xyz(1)-10e-6_dp)<tol) then
        ir=ir+1;right(ir)=i
      end if
      if(abs(xyz(2))<tol) then
        isym=isym+1;symmetry(isym)=i
      end if
    end do
    allocate(bc(11))
    bc(1)%field_id=F_U1;bc(1)%node_ids=[a]
    bc(2)%field_id=F_U2;bc(2)%node_ids=symmetry
    bc(3)%field_id=F_U3;bc(3)%node_ids=[a,b]
    do i=4,7
      bc(2*i-4)%field_id=i;bc(2*i-4)%node_ids=left
      bc(2*i-3)%field_id=i;bc(2*i-3)%node_ids=right
    end do
    bc(5)%value_end=voltage
  end subroutine define_half_y_boundaries

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

  subroutine write_nodes(directory,mesh,material,state,status,filename)
    character(len=*), intent(in) :: directory
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    type(SolverStatusType), intent(out) :: status
    character(len=*), intent(in), optional :: filename
    character(len=80) :: name
    integer :: unit_id,io,node
    real(dp) :: nratio,pratio
    call clear_status(status)
    name='nodes.csv'
    if(present(filename)) name=filename
    open(newunit=unit_id,file=trim(directory)//'/'//trim(name),status='new',action='write',iostat=io)
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

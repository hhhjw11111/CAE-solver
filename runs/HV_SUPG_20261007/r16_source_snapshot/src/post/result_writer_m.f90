module result_writer_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : F_U1, F_U2, F_U3, F_PHI, F_ETA_N, F_ETA_P, &
    F_THETA, global_dof
  use build_info_m, only : SOLVER_NAME, SOLVER_VERSION, THEORY_VERSION, ARCHITECTURE_ID
  use input_schema_m, only : CaseConfigType
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INPUT_ERROR
  use continuation_solver_m, only : ContinuationReportType, CONTINUATION_PHASE_START, &
    CONTINUATION_PHASE_VOLTAGE, CONTINUATION_PHASE_JOULE
  use independent_flux_audit_m, only : IndependentFluxReportType
  implicit none
  private

  public :: write_solution_result

contains

  subroutine write_solution_result(result_path,config,mesh,material,state,load,report, &
      flux_audit,flux_audit4,status)
    character(len=*), intent(in) :: result_path
    type(CaseConfigType), intent(in) :: config
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(ContinuationReportType), intent(in) :: report
    type(IndependentFluxReportType), intent(in) :: flux_audit,flux_audit4
    type(SolverStatusType), intent(out) :: status
    integer :: unit_id,io_status,point,node,component
    character(len=512) :: io_message
    real(dp) :: values(12),eta_n,eta_p,n,p
    real(dp) :: volume_scale,power_scale,quadrature_volume_error
    real(dp) :: quadrature_joule_error,quadrature_heat_error

    call clear_status(status)
    if (len_trim(result_path) == 0) then
      call set_status(status,STATUS_INPUT_ERROR,'Result path is empty.')
      return
    end if
    open(newunit=unit_id,file=trim(result_path),status='replace',action='write', &
      iostat=io_status,iomsg=io_message)
    if (io_status /= 0) then
      call set_status(status,STATUS_INPUT_ERROR,'Cannot open result file: '//trim(io_message))
      return
    end if

    write(unit_id,'(a)') 'solver='//SOLVER_NAME
    write(unit_id,'(a)') 'version='//SOLVER_VERSION
    write(unit_id,'(a)') 'theory_version='//THEORY_VERSION
    write(unit_id,'(a)') 'architecture_id='//ARCHITECTURE_ID
    write(unit_id,'(a)') 'run_result=PASS'
    write(unit_id,'(a)') 'run_result_scope=nonlinear_solve_and_output_not_flux_acceptance'
    write(unit_id,'(a)') 'case_id='//trim(config%case_id)
    write(unit_id,'(a)') 'model_level='//trim(config%model_level)
    write(unit_id,'(a)') 'problem_type='//trim(config%problem_type)
    write(unit_id,'(a)') 'material_profile='//trim(config%material_profile)
    write(unit_id,'(a,es24.16e3)') 'hole_radius=',config%hole_radius
    write(unit_id,'(a,i0)') 'node_count=',mesh%number_of_nodes()
    write(unit_id,'(a,i0)') 'element_count=',mesh%number_of_elements()
    write(unit_id,'(a,i0)') 'accepted_points=',report%accepted_points
    write(unit_id,'(a,i0)') 'rejected_attempts=',report%rejected_attempts
    write(unit_id,'(a,es24.16e3)') 'final_voltage_factor=',load%voltage_factor
    write(unit_id,'(a,es24.16e3)') 'final_joule_factor=',load%joule_factor
    write(unit_id,'(a,es24.16e3)') 'maximum_balance_defect=',report%maximum_balance_defect
    write(unit_id,'(a,es24.16e3)') 'maximum_boundary_error=',report%maximum_boundary_error
    volume_scale=max(abs(flux_audit%volume),1.0e-300_dp)
    power_scale=max(abs(flux_audit%joule),abs(flux_audit%total_current_work), &
      abs(flux_audit%terminal_power),sum(abs(flux_audit%heat_out)),tiny(1.0_dp))
    quadrature_volume_error=abs(flux_audit%volume-flux_audit4%volume)/volume_scale
    quadrature_joule_error=abs(flux_audit%joule-flux_audit4%joule)/power_scale
    quadrature_heat_error=maxval(abs(flux_audit%heat_out-flux_audit4%heat_out))/power_scale
    write(unit_id,'(a)') 'independent_flux_audit=COMPUTED'
    write(unit_id,'(a,i0)') 'independent_flux_audit_order=',5
    write(unit_id,'(a,i0)') 'independent_flux_exterior_faces=',flux_audit%exterior_faces
    write(unit_id,'(a,i0)') 'independent_flux_internal_faces=',flux_audit%internal_faces
    write(unit_id,'(a)') 'independent_flux_status=COMPUTED_NOT_ACCEPTANCE_TESTED'
    write(unit_id,'(a)') 'independent_flux_face_evaluation=EXACT_BOUNDARY_NO_OFFSET'
    write(unit_id,'(a)') 'independent_flux_sets=min_x_contact,max_x_contact,other_exterior'
    write(unit_id,'(a)') 'independent_flux_shared_kernels=shape_geometry_material_point'
    write(unit_id,'(a,es24.16e3)') 'independent_flux_power_normalization_W=',power_scale
    write(unit_id,'(a,es24.16e3)') 'independent_flux_volume=',flux_audit%volume
    write(unit_id,'(a,es24.16e3)') 'independent_flux_area=',flux_audit%area
    write(unit_id,'(a,3(es24.16e3,a))') 'independent_flux_area_normal=', &
      flux_audit%area_normal(1),',',flux_audit%area_normal(2),',',flux_audit%area_normal(3),','
    write(unit_id,'(a,3(es24.16e3,a))') 'independent_electron_out_A=', &
      flux_audit%electron_out(1),',',flux_audit%electron_out(2),',',flux_audit%electron_out(3),','
    write(unit_id,'(a,3(es24.16e3,a))') 'independent_hole_out_A=', &
      flux_audit%hole_out(1),',',flux_audit%hole_out(2),',',flux_audit%hole_out(3),','
    write(unit_id,'(a,3(es24.16e3,a))') 'independent_heat_out_W=', &
      flux_audit%heat_out(1),',',flux_audit%heat_out(2),',',flux_audit%heat_out(3),','
    write(unit_id,'(a,es24.16e3)') 'independent_electric_displacement_out_C=', &
      flux_audit%electric_displacement_out
    write(unit_id,'(a,es24.16e3)') 'independent_charge_C=',flux_audit%charge
    write(unit_id,'(a,es24.16e3)') 'independent_recombination_charge_A=', &
      flux_audit%recombination_charge
    write(unit_id,'(a,es24.16e3)') 'independent_joule_W=',flux_audit%joule
    write(unit_id,'(a,es24.16e3)') 'independent_terminal_power_W=',flux_audit%terminal_power
    write(unit_id,'(a,es24.16e3)') 'independent_total_current_work_W=', &
      flux_audit%total_current_work
    write(unit_id,'(a,es24.16e3)') 'independent_diffusion_work_W=',flux_audit%diffusion_work
    write(unit_id,'(a,es24.16e3)') 'independent_absolute_diffusion_work_W=', &
      flux_audit%absolute_diffusion_work
    write(unit_id,'(a,es24.16e3)') 'independent_heat_balance_defect_W=', &
      sum(flux_audit%heat_out)-flux_audit%joule
    write(unit_id,'(a,es24.16e3)') 'independent_heat_balance_relative=', &
      abs(sum(flux_audit%heat_out)-flux_audit%joule)/power_scale
    write(unit_id,'(a,es24.16e3)') 'independent_charge_balance_C=', &
      flux_audit%electric_displacement_out-flux_audit%charge
    write(unit_id,'(a,es24.16e3)') 'independent_power_balance_relative=', &
      abs(flux_audit%terminal_power-flux_audit%total_current_work)/power_scale
    write(unit_id,'(a,es24.16e3)') 'independent_quadrature_volume_error=', &
      quadrature_volume_error
    write(unit_id,'(a,es24.16e3)') 'independent_quadrature_joule_error=', &
      quadrature_joule_error
    write(unit_id,'(a,es24.16e3)') 'independent_quadrature_heat_error=', &
      quadrature_heat_error
    write(unit_id,'(a)') 'continuation_columns=point,phase,parameter,voltage_factor,'// &
      'joule_factor,newton_iterations,residual_rms,update_rms,charge_balance,'// &
      'electron_balance,hole_balance,heat_balance'
    do point=1,report%accepted_points
      write(unit_id,'(i0,a,a,a,3(es24.16e3,a),i0,a,6(es24.16e3,a))') point,',', &
        trim(phase_name(report%phase(point))),',',report%path_parameter(point),',', &
        report%voltage_factor(point),',',report%joule_factor(point),',', &
        report%newton_iterations(point),',',report%residual_rms(point),',', &
        report%update_rms(point),',',report%balance(1,point),',',report%balance(2,point),',', &
        report%balance(3,point),',',report%balance(4,point),','
    end do
    write(unit_id,'(a)') 'node_columns=node_id,x,y,z,u1,u2,u3,phi,eta_n,eta_p,theta,n,p'
    do node=1,mesh%number_of_nodes()
      eta_n=state%x(global_dof(node,F_ETA_N))
      eta_p=state%x(global_dof(node,F_ETA_P))
      n=material%reference_n*exp(eta_n)
      p=material%reference_p*exp(eta_p)
      values(1:3)=mesh%coordinates(:,node)
      values(4)=state%x(global_dof(node,F_U1))
      values(5)=state%x(global_dof(node,F_U2))
      values(6)=state%x(global_dof(node,F_U3))
      values(7)=state%x(global_dof(node,F_PHI))
      values(8)=eta_n
      values(9)=eta_p
      values(10)=state%x(global_dof(node,F_THETA))
      values(11)=n
      values(12)=p
      if (.not. all(ieee_is_finite(values)) .or. n <= 0.0_dp .or. p <= 0.0_dp) then
        close(unit_id)
        call set_status(status,STATUS_INPUT_ERROR, &
          'Result contains a non-finite or non-positive nodal value.')
        return
      end if
      write(unit_id,'(i0,12(a,es24.16e3))') node, &
        (',',values(component),component=1,12)
    end do
    close(unit_id,iostat=io_status,iomsg=io_message)
    if (io_status /= 0) then
      call set_status(status,STATUS_INPUT_ERROR,'Cannot close result file: '//trim(io_message))
    end if
  end subroutine write_solution_result

  pure function phase_name(phase) result(name)
    integer, intent(in) :: phase
    character(len=8) :: name

    select case(phase)
    case(CONTINUATION_PHASE_START)
      name='START'
    case(CONTINUATION_PHASE_VOLTAGE)
      name='VOLTAGE'
    case(CONTINUATION_PHASE_JOULE)
      name='JOULE'
    case default
      name='UNKNOWN'
    end select
  end function phase_name

end module result_writer_m

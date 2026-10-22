module production_runner_m
  use, intrinsic :: iso_fortran_env, only : output_unit
  use kinds_m, only : dp
  use input_schema_m, only : CaseConfigType, MODEL_M2
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INPUT_ERROR
  use native_input_m, only : read_native_case
  use mesh_checks_m, only : validate_mesh
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context
  use dirichlet_m, only : DirichletBCType
  use reference_bar_problem_m, only : initialize_voltage_driven_problem
  use continuation_solver_m, only : ContinuationReportType, solve_load_continuation
  use independent_flux_audit_m, only : IndependentFluxReportType, &
    integrate_independent_fluxes
  use result_writer_m, only : write_solution_result
  implicit none
  private

  public :: run_production_case

contains

  subroutine run_production_case(case_path,result_path,status)
    character(len=*), intent(in) :: case_path,result_path
    type(SolverStatusType), intent(out) :: status
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    type(MaterialType), allocatable :: materials(:)
    type(StateType) :: state
    type(LoadStateType) :: current_load,target_load
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(SparseAssemblyContextType) :: context
    type(ContinuationReportType) :: report
    type(IndependentFluxReportType) :: flux_audit,flux_audit4
    type(SolverStatusType) :: local_status
    real(dp) :: minimum_determinant,total_volume,volume_scale,power_scale
    real(dp) :: quadrature_volume_error,quadrature_joule_error
    real(dp) :: quadrature_heat_error

    call clear_status(status)
    call read_native_case(case_path,config,mesh,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    if (trim(config%model_level) /= MODEL_M2) then
      call set_status(status,STATUS_INPUT_ERROR, &
        'Production run currently supports only the frozen M2 model level.')
      return
    end if
    call validate_mesh(mesh,minimum_determinant,total_volume,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call initialize_voltage_driven_problem(config,mesh,materials,state,current_load, &
      boundary_conditions,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if

    target_load%voltage_factor=config%target_voltage_factor
    target_load%joule_factor=config%target_joule_factor
    call solve_load_continuation(mesh,materials,context,boundary_conditions, &
      config%field_scales,config%nonlinear_options,config%continuation_options, &
      target_load,state,current_load,report,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if
    call integrate_independent_fluxes(mesh,materials,state,current_load,5,flux_audit,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,local_status%code, &
        'Independent order-5 flux audit failed: '//trim(local_status%message)); return
    end if
    call integrate_independent_fluxes(mesh,materials,state,current_load,4,flux_audit4,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,local_status%code, &
        'Independent order-4 flux audit failed: '//trim(local_status%message)); return
    end if
    volume_scale=max(total_volume,1.0e-300_dp)
    power_scale=max(abs(flux_audit%joule),abs(flux_audit%total_current_work), &
      abs(flux_audit%terminal_power),sum(abs(flux_audit%heat_out)),tiny(1.0_dp))
    quadrature_volume_error=abs(flux_audit%volume-flux_audit4%volume)/volume_scale
    quadrature_joule_error=abs(flux_audit%joule-flux_audit4%joule)/power_scale
    quadrature_heat_error=maxval(abs(flux_audit%heat_out-flux_audit4%heat_out))/ &
      power_scale
    call write_solution_result(result_path,config,mesh,materials(1),state,current_load, &
      report,flux_audit,flux_audit4,local_status)
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status); return
    end if

    write(output_unit,'(a)') 'run_result=PASS'
    write(output_unit,'(a)') 'run_result_scope=nonlinear_solve_and_output_not_flux_acceptance'
    write(output_unit,'(a)') 'case_id='//trim(config%case_id)
    write(output_unit,'(a,i0)') 'node_count=',mesh%number_of_nodes()
    write(output_unit,'(a,i0)') 'element_count=',mesh%number_of_elements()
    write(output_unit,'(a,es24.16e3)') 'minimum_jacobian_determinant=',minimum_determinant
    write(output_unit,'(a,es24.16e3)') 'integrated_volume=',total_volume
    write(output_unit,'(a,i0)') 'accepted_points=',report%accepted_points
    write(output_unit,'(a,i0)') 'rejected_attempts=',report%rejected_attempts
    write(output_unit,'(a,i0)') 'total_newton_iterations=',report%total_newton_iterations
    write(output_unit,'(a,es24.16e3)') 'maximum_balance_defect=', &
      report%maximum_balance_defect
    write(output_unit,'(a,es24.16e3)') 'maximum_boundary_error=', &
      report%maximum_boundary_error
    write(output_unit,'(a,i0)') 'independent_flux_audit_order=',5
    write(output_unit,'(a,i0)') 'independent_flux_exterior_faces=',flux_audit%exterior_faces
    write(output_unit,'(a,i0)') 'independent_flux_internal_faces=',flux_audit%internal_faces
    write(output_unit,'(a)') 'independent_flux_status=COMPUTED_NOT_ACCEPTANCE_TESTED'
    write(output_unit,'(a)') 'independent_flux_face_evaluation=EXACT_BOUNDARY_NO_OFFSET'
    write(output_unit,'(a,es24.16e3)') 'independent_flux_power_normalization_W=',power_scale
    write(output_unit,'(a,es24.16e3)') 'independent_flux_joule_W=',flux_audit%joule
    write(output_unit,'(a,es24.16e3)') 'independent_flux_terminal_power_W=', &
      flux_audit%terminal_power
    write(output_unit,'(a,es24.16e3)') 'independent_flux_heat_balance_defect_W=', &
      sum(flux_audit%heat_out)-flux_audit%joule
    write(output_unit,'(a,es24.16e3)') 'independent_flux_heat_balance_relative=', &
      abs(sum(flux_audit%heat_out)-flux_audit%joule)/power_scale
    write(output_unit,'(a,es24.16e3)') 'independent_flux_power_balance_relative=', &
      abs(flux_audit%terminal_power-flux_audit%total_current_work)/power_scale
    write(output_unit,'(a,es24.16e3)') 'independent_flux_charge_balance_C=', &
      flux_audit%electric_displacement_out-flux_audit%charge
    write(output_unit,'(a,es24.16e3)') 'independent_flux_quadrature_volume_error=', &
      quadrature_volume_error
    write(output_unit,'(a,es24.16e3)') 'independent_flux_quadrature_joule_error=', &
      quadrature_joule_error
    write(output_unit,'(a,es24.16e3)') 'independent_flux_quadrature_heat_error=', &
      quadrature_heat_error
    write(output_unit,'(a,es24.16e3)') 'final_voltage_factor=',current_load%voltage_factor
    write(output_unit,'(a,es24.16e3)') 'final_joule_factor=',current_load%joule_factor
  end subroutine run_production_case

  subroutine copy_failure(source,destination)
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,source%code,trim(source%message))
  end subroutine copy_failure

end module production_runner_m

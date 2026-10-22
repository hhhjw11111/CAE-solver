module temperature_control_check_m
  use kinds_m, only: dp
  use field_layout_m, only: NDOF_NODE,F_PHI,F_THETA,global_dof
  use input_schema_m, only: CaseConfigType,PROBLEM_REFERENCE_BAR,MATERIAL_TPS_G08_REFERENCE
  use mesh_types_m, only: MeshType
  use model_types_m, only: MaterialType,StateType,LoadStateType,DiagnosticsType
  use status_m, only: SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  use sparse_pattern_m, only: SparseAssemblyContextType,initialize_sparse_context
  use dirichlet_m, only: DirichletBCType,build_fixed_dof_mask
  use reference_bar_problem_m, only: initialize_reference_bar_problem
  use global_assembly_check_m, only: build_two_element_mesh
  use scaling_m, only: build_variable_scaling
  use continuation_solver_m, only: ContinuationReportType,solve_load_continuation
  use newton_solver_m, only: NewtonReportType,solve_single_load_newton
  use temperature_control_m, only: solve_temperature_controlled_newton
  use balance_audit_m, only: audit_global_balances
  implicit none
  private
  public :: verify_temperature_control
contains
  subroutine verify_temperature_control(status)
    type(SolverStatusType), intent(out) :: status
    type(SolverStatusType) :: local_status
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    type(MaterialType), allocatable :: materials(:)
    type(StateType) :: reference,inverse,replay
    type(LoadStateType) :: load,target,inverse_load
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: bc(:)
    type(ContinuationReportType) :: continuation
    type(NewtonReportType) :: report
    type(DiagnosticsType) :: diagnostics
    real(dp), allocatable :: scales(:)
    logical, allocatable :: fixed(:)
    real(dp) :: theta,field_error,load_error,balance(4),bc_error
    integer :: pass,node,control,equation
    call clear_status(status)
    do pass=1,2
      config=CaseConfigType()
      config%problem_type=PROBLEM_REFERENCE_BAR;config%material_profile=MATERIAL_TPS_G08_REFERENCE
      config%target_voltage=-10
      config%field_scales=[1e-10_dp,1e-10_dp,1e-10_dp,10.0_dp,.1_dp,.1_dp,10.0_dp]
      config%nonlinear_options%max_newton=30
      config%nonlinear_options%residual_tolerance=1e-9_dp;config%nonlinear_options%update_tolerance=1e-9_dp
      config%continuation_options%maximum_newton_per_attempt=30
      call build_two_element_mesh(mesh);mesh%coordinates=1e-6_dp*mesh%coordinates
      call initialize_reference_bar_problem(config,mesh,materials,reference,load,bc,status)
      if(.not.status%is_ok()) return
      materials(1)%transport_supg=pass==2;materials(1)%supg_directional_tau=.true.
      call initialize_sparse_context(mesh,context,status)
      if(.not.status%is_ok()) return
      target%voltage_factor=.8_dp;target%joule_factor=1
      call solve_load_continuation(mesh,materials,context,bc,config%field_scales,config%nonlinear_options, &
        config%continuation_options,target,reference,load,continuation,status)
      if(.not.status%is_ok()) return
      call build_fixed_dof_mask(bc,mesh%number_of_nodes(),fixed,status)
      if(.not.status%is_ok()) return
      theta=-huge(1.0_dp);control=0
      do node=1,mesh%number_of_nodes()
        equation=global_dof(node,F_THETA)
        if(fixed(equation)) cycle
        ! Interior of the first element: the opposite electrode is outside its
        ! original CSR neighbors, so the new parameter-column insertion is tested.
        if(mesh%coordinates(1,node)<.4e-6_dp.or.mesh%coordinates(1,node)>.6e-6_dp) cycle
        if(reference%x(equation)>theta) then
          theta=reference%x(equation);control=node
        end if
      end do
      inverse=reference;inverse_load=load;inverse_load%voltage_factor=.75_dp
      inverse%x(F_PHI::NDOF_NODE)=inverse%x(F_PHI::NDOF_NODE)*(.75_dp/.8_dp)
      inverse%x(F_THETA::NDOF_NODE)=inverse%x(F_THETA::NDOF_NODE)*(.75_dp/.8_dp)**2
      call solve_temperature_controlled_newton(mesh,materials,inverse_load,context,bc,config%field_scales, &
        config%nonlinear_options,control,theta,inverse,report,status)
      if(.not.status%is_ok()) return
      call build_variable_scaling(mesh%number_of_nodes(),config%field_scales,scales,status)
      if(.not.status%is_ok()) return
      field_error=maxval(abs((inverse%x-reference%x)/scales))
      load_error=abs(inverse_load%voltage_factor-target%voltage_factor)
      call audit_global_balances(mesh,materials,inverse,inverse_load,context,bc,balance,bc_error,diagnostics,status)
      if(.not.status%is_ok()) return
      write(*,'(a,i0,a,3(es24.16e3,a))') 'temperature_control_test=',pass, &
        ',load_error=',load_error,',field_error=',field_error,',balance=',maxval(balance),','
      if(load_error>=1e-9_dp.or.field_error>=1e-8_dp.or.maxval(balance)>=1e-8_dp.or.bc_error>=1e-10_dp) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Inverse parameterization disagrees with the ordinary voltage solution.')
        return
      end if
      replay=inverse
      call solve_single_load_newton(mesh,materials,inverse_load,context,bc,config%field_scales, &
        config%nonlinear_options,replay,report,status)
      if(.not.status%is_ok()) return
      if(maxval(abs((replay%x-inverse%x)/scales))>=1e-8_dp) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Fixed-voltage replay changed the inverse state.');return
      end if
      call solve_temperature_controlled_newton(mesh,materials,inverse_load,context,bc,config%field_scales, &
        config%nonlinear_options,bc(10)%node_ids(1),theta,inverse,report,local_status)
      if(local_status%is_ok()) then
        call set_status(status,STATUS_VERIFICATION_FAILED,'Fixed temperature node was incorrectly accepted as the control.');return
      end if
    end do
    write(*,'(a)') 'verification_id=TCTRL_original_PDE_and_voltage_replay'
    write(*,'(a)') 'verification_result=PASS'
  end subroutine verify_temperature_control
end module temperature_control_check_m

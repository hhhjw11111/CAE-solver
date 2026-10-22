module arclength_solver_check_m
  use kinds_m, only : dp
  use model_types_m, only : StateType,LoadStateType,NonlinearOptionsType,ContinuationOptionsType
  use sparse_pattern_m, only : SparseAssemblyContextType
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  use arclength_solver_m, only : trace_arclength_system
  use field_layout_m, only : NDOF_NODE
  use mesh_types_m, only : MeshType
  use continuation_solver_m, only : ContinuationReportType
  implicit none
  private
  public :: verify_arclength_saddle_node
contains
  subroutine verify_arclength_saddle_node(status,state_unit)
    type(SolverStatusType), intent(out) :: status
    real(dp), intent(in), optional :: state_unit
    type(SparseAssemblyContextType) :: context
    type(StateType) :: state
    type(LoadStateType) :: load
    type(ContinuationReportType) :: report
    type(NonlinearOptionsType) :: options
    type(ContinuationOptionsType) :: arc_options
    real(dp) :: scales(1)
    real(dp) :: unit_scale
    logical :: fixed(1),limited(1)
    call clear_status(status)
    unit_scale=1.0_dp
    if(present(state_unit)) unit_scale=state_unit
    context%number_of_equations=1;context%number_of_nonzeros=1
    allocate(context%ia(2),context%ja(1));context%ia=[1,2];context%ja=[1]
    allocate(state%x(1));state%x=-unit_scale;load%voltage_factor=0.0_dp;load%joule_factor=0.0_dp
    scales=unit_scale;fixed=.false.;limited=.false.
    options%max_newton=20;options%residual_tolerance=1.0e-10_dp;options%update_tolerance=1.0e-10_dp
    options%armijo_c1=1.0e-4_dp;options%backtrack_factor=0.5_dp;options%minimum_step=1.0e-8_dp
    options%maximum_log_density_update=10.0_dp
    arc_options=ContinuationOptionsType();arc_options%arclength_enabled=.true.
    arc_options%initial_step=0.2_dp;arc_options%minimum_step=1.0e-4_dp;arc_options%maximum_step=0.25_dp
    arc_options%maximum_attempts=80;arc_options%maximum_newton_per_attempt=15
    arc_options%arclength_psi_state=1.0_dp;arc_options%arclength_psi_voltage=1.0_dp
    arc_options%arclength_fold_tolerance=2.0e-3_dp;arc_options%arclength_maximum_refinements=20
    call trace_arclength_system(context,scales,fixed,limited,options,arc_options,1.5_dp, &
      saddle_node_residual,state,load,report,status)
    write(*,'(a)') 'verification_id=ARC_saddle_node'
    write(*,'(a,l1)') 'arc_fold_detected=',report%arclength_fold_detected
    write(*,'(a,l1)') 'arc_fold_located=',report%arclength_fold_located
    write(*,'(a,es24.16e3)') 'arc_fold_voltage_factor=',report%arclength_fold_voltage_factor
    write(*,'(a,es24.16e3)') 'arc_fold_bracket_width=',report%arclength_fold_bracket_width
    write(*,'(a,es24.16e3)') 'arc_final_state=',state%x(1)
    if(present(state_unit)) write(*,'(a,es24.16e3)') 'arc_state_unit=',unit_scale
    if(.not.status%is_ok()) return
    if(.not.report%arclength_fold_detected.or..not.report%arclength_fold_located.or. &
        abs(report%arclength_fold_voltage_factor-1.0_dp)>5.0e-3_dp.or. &
        report%arclength_fold_bracket_width>arc_options%arclength_fold_tolerance) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'Pseudo-arclength saddle-node location failed.')
      write(*,'(a)') 'verification_result=FAIL';return
    end if
    write(*,'(a)') 'verification_result=PASS'
  contains
    subroutine saddle_node_residual(x,lambda,need_jacobian,residual,jacobian,parameter,status)
      real(dp), intent(in) :: x(:),lambda
      logical, intent(in) :: need_jacobian
      real(dp), intent(out) :: residual(:),jacobian(:),parameter(:)
      type(SolverStatusType), intent(out) :: status
      call clear_status(status)
      residual(1)=(x(1)/unit_scale)**2+lambda-1.0_dp
      jacobian(1)=2.0_dp*x(1)/unit_scale**2;parameter(1)=1.0_dp
      if(.not.need_jacobian) jacobian=0.0_dp
    end subroutine
  end subroutine verify_arclength_saddle_node
end module arclength_solver_check_m

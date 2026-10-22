module balance_audit_m
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE, F_PHI, F_ETA_N, F_ETA_P, F_THETA
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, DiagnosticsType
  use status_m, only : SolverStatusType, clear_status
  use sparse_pattern_m, only : SparseAssemblyContextType
  use dirichlet_m, only : DirichletBCType, project_dirichlet_state, build_fixed_dof_mask
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, &
    assemble_global_system
  implicit none
  private

  public :: audit_global_balances

contains

  subroutine audit_global_balances(mesh,materials,state,load,context,boundary_conditions, &
                                   balance,boundary_error,diagnostics,status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    real(dp), intent(out) :: balance(4)
    real(dp), intent(out) :: boundary_error
    type(DiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    type(StateType) :: projected_state
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: unconstrained_system
    logical, allocatable :: fixed(:)
    integer :: fields(4),group

    call clear_status(status)
    request%need_jacobian=.false.
    call assemble_global_system(mesh,materials,state,load,context,request=request, &
      system=unconstrained_system,diagnostics=diagnostics,status=status)
    if (.not. status%is_ok()) return
    call build_fixed_dof_mask(boundary_conditions,mesh%number_of_nodes(),fixed,status)
    if (.not. status%is_ok()) return
    fields=[F_PHI,F_ETA_N,F_ETA_P,F_THETA]
    do group=1,4
      balance(group)=field_free_balance(unconstrained_system%residual,fixed,fields(group))
    end do
    allocate(projected_state%x(size(state%x)))
    projected_state%x=state%x
    call project_dirichlet_state(projected_state,boundary_conditions,load, &
      mesh%number_of_nodes(),status)
    if (.not. status%is_ok()) return
    boundary_error=maxval(abs(projected_state%x-state%x))
  end subroutine audit_global_balances

  real(dp) function field_free_balance(residual,fixed,field) result(value)
    real(dp), intent(in) :: residual(:)
    logical, intent(in) :: fixed(:)
    integer, intent(in) :: field
    real(dp) :: free_sum,absolute_sum
    integer :: equation

    free_sum=0.0_dp
    absolute_sum=0.0_dp
    do equation=field,size(residual),NDOF_NODE
      absolute_sum=absolute_sum+abs(residual(equation))
      if (.not. fixed(equation)) free_sum=free_sum+residual(equation)
    end do
    value=abs(free_sum)/max(absolute_sum,tiny(1.0_dp))
  end function field_free_balance

end module balance_audit_m

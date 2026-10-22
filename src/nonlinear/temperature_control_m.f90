! Diagnostic reparameterization: prescribe one free nodal temperature and solve
! for the voltage load factor. Every original PDE row is retained. The prescribed
! temperature's column is replaced by the derivative of the electrode constraint.
! Returned points require an ordinary fixed-voltage replay and independent flux audit.
module temperature_control_m
  use, intrinsic :: iso_fortran_env, only: output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use kinds_m, only: dp
  use field_layout_m, only: NDOF_NODE,F_PHI,F_THETA,F_ETA_N,F_ETA_P,global_dof
  use mesh_types_m, only: MeshType
  use model_types_m, only: MaterialType,StateType,LoadStateType,NonlinearOptionsType,DiagnosticsType
  use status_m, only: SolverStatusType,clear_status,set_status,STATUS_INVALID_ARGUMENT,STATUS_NONLINEAR_SOLVER_FAILED
  use sparse_pattern_m, only: SparseAssemblyContextType,csr_matvec
  use dirichlet_m, only: DirichletBCType,project_dirichlet_state,build_fixed_dof_mask
  use global_assembly_m, only: GlobalAssemblyRequestType,GlobalSystemType,assemble_global_system
  use scaling_m, only: build_variable_scaling,equilibrate_csr_system
  use pardiso_solver_m, only: PardisoDiagnosticsType,solve_csr_pardiso
  use newton_solver_m, only: NewtonReportType
  implicit none
  private
  public :: solve_temperature_controlled_newton
contains
  subroutine solve_temperature_controlled_newton(mesh,materials,load,base,bc,field_scales,options, &
                                                 control_node,theta_goal,state,report,status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(LoadStateType), intent(inout) :: load
    type(SparseAssemblyContextType), intent(in) :: base
    type(DirichletBCType), intent(in) :: bc(:)
    real(dp), intent(in) :: field_scales(NDOF_NODE),theta_goal
    type(NonlinearOptionsType), intent(in) :: options
    integer, intent(in) :: control_node
    type(StateType), intent(inout) :: state
    type(NewtonReportType), intent(out) :: report
    type(SolverStatusType), intent(out) :: status
    type(SparseAssemblyContextType) :: ctx
    type(StateType) :: trial
    type(LoadStateType) :: trial_load
    type(GlobalAssemblyRequestType) :: request,trial_request
    type(GlobalSystemType) :: system,trial_system
    type(DiagnosticsType) :: diagnostics,trial_diagnostics
    type(PardisoDiagnosticsType) :: linear_diagnostics
    type(SolverStatusType) :: local_status
    logical, allocatable :: fixed(:)
    integer, allocatable :: mapping(:),control_positions(:)
    real(dp), allocatable :: parameter(:),scales(:),values(:),sv(:),srhs(:),rs(:),dx(:),delta(:),action(:),history(:)
    real(dp), allocatable :: rhs(:),trial_scaled(:)
    real(dp) :: residual,update,alpha,limit,logmax,slope,trial_residual
    integer :: n,control,row,b,j,iteration,p
    logical :: accepted
    call clear_status(status);report=NewtonReportType()
    n=base%number_of_equations
    if(control_node<1.or.control_node>mesh%number_of_nodes().or..not.ieee_is_finite(theta_goal).or. &
       theta_goal<=0.or.load%joule_factor/=1.0_dp.or.options%max_newton<1.or. &
       options%minimum_step<=0.or.options%backtrack_factor<=0.or.options%backtrack_factor>=1.or. &
       options%maximum_log_density_update<=0.or.options%residual_tolerance<=0.or.options%update_tolerance<=0.or. &
       options%armijo_c1<=0.or.options%armijo_c1>=.5_dp.or..not.ieee_is_finite(load%voltage_factor)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Temperature control requires a positive free temperature and full Joule heat.')
      return
    end if
    if(options%trust_region) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Trust region is supported only by the fixed-load Newton solver.');return
    end if
    control=global_dof(control_node,F_THETA)
    call build_fixed_dof_mask(bc,mesh%number_of_nodes(),fixed,status)
    if(.not.status%is_ok()) return
    if(fixed(control)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Temperature control node is already a Dirichlet temperature.');return
    end if
    allocate(parameter(n));parameter=0
    do b=1,size(bc)
      if(bc(b)%field_id/=F_PHI.and.bc(b)%value_end/=bc(b)%value_start) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'Temperature control only supports voltage-dependent electrical contacts.')
        return
      end if
      if(bc(b)%field_id/=F_PHI) cycle
      do j=1,size(bc(b)%node_ids)
        row=global_dof(bc(b)%node_ids(j),F_PHI)
        parameter(row)=-(bc(b)%value_end-bc(b)%value_start)
      end do
    end do
    if(maxval(abs(parameter))==0) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Temperature control requires a nonzero voltage boundary range.');return
    end if
    call parameter_context(base,control,parameter,ctx,mapping,control_positions)
    call build_variable_scaling(mesh%number_of_nodes(),field_scales,scales,status)
    if(.not.status%is_ok()) return
    scales(control)=1.0_dp ! This linear-system coordinate is lambda, not theta.
    allocate(values(ctx%number_of_nonzeros),dx(n),delta(n),action(n),trial%x(n),history(options%max_newton+1), &
      rhs(n),trial_scaled(n))
    history=huge(1.0_dp);update=huge(1.0_dp)
    state%x(control)=theta_goal
    call project_dirichlet_state(state,bc,load,mesh%number_of_nodes(),status)
    if(.not.status%is_ok()) return
    request%need_jacobian=.true.;trial_request%need_jacobian=.false.
    do iteration=0,options%max_newton
      call assemble_global_system(mesh,materials,state,load,base,bc,request,system,diagnostics,status)
      if(.not.status%is_ok()) return
      values=0
      do p=1,base%number_of_nonzeros
        values(mapping(p))=system%jacobian_values(p)
      end do
      do row=1,n
        if(control_positions(row)>0) values(control_positions(row))=parameter(row)
      end do
      rhs=-system%residual
      call equilibrate_csr_system(ctx,values,rhs,scales,sv,srhs,rs,status)
      if(.not.status%is_ok()) return
      residual=free_rms(srhs,fixed);history(iteration+1)=residual
      report%iterations=iteration;report%final_residual_rms=residual;report%final_diagnostics=diagnostics
      if(iteration==0) report%initial_residual_rms=residual
      if(options%verbose) then
        write(output_unit,'(a,i0,a,3(es16.8,a))') 'temperature_newton_iteration=',iteration, &
          ',voltage_factor=',load%voltage_factor,',theta_control=',theta_goal,',residual=',residual,','
        flush(output_unit)
      end if
      if(residual<options%residual_tolerance.and.(iteration==0.or.update<options%update_tolerance)) then
        report%converged=.true.;report%final_update_rms=update
        if(iteration==0) report%final_update_rms=0
        return
      end if
      if(iteration==options%max_newton) exit
      if(options%stagnation_window>0.and.residual>options%residual_tolerance) then
        if(iteration>=options%stagnation_window) then
          if(residual>=.99_dp*history(iteration+1-options%stagnation_window)) exit
        end if
      end if
      call solve_csr_pardiso(n,ctx%ia,ctx%ja,sv,srhs,dx,linear_diagnostics,status)
      if(.not.status%is_ok()) return
      report%total_linear_solves=report%total_linear_solves+1
      delta=dx*scales
      logmax=max(maxval(abs(delta(F_ETA_N::NDOF_NODE))),maxval(abs(delta(F_ETA_P::NDOF_NODE))))
      limit=min(1.0_dp,options%maximum_log_density_update/max(logmax,tiny(1.0_dp)))
      dx=limit*dx;delta=limit*delta
      call csr_matvec(ctx,sv,dx,action,status)
      if(.not.status%is_ok()) return
      slope=0
      do row=1,n
        if(.not.fixed(row)) slope=slope-srhs(row)*action(row)
      end do
      slope=slope/real(count(.not.fixed),dp)
      if(.not.ieee_is_finite(slope).or.slope>=0) exit
      alpha=1;accepted=.false.
      do while(alpha>=options%minimum_step)
        trial_load=load;trial_load%voltage_factor=load%voltage_factor+alpha*delta(control)
        if(ieee_is_finite(trial_load%voltage_factor).and.trial_load%voltage_factor>0.and. &
           trial_load%voltage_factor<1.25_dp) then
          trial%x=state%x+alpha*delta;trial%x(control)=theta_goal
          call project_dirichlet_state(trial,bc,trial_load,mesh%number_of_nodes(),local_status)
          if(local_status%is_ok()) then
            call assemble_global_system(mesh,materials,trial,trial_load,base,bc,trial_request, &
              trial_system,trial_diagnostics,local_status)
            if(local_status%is_ok()) then
              trial_scaled=trial_system%residual/rs
              trial_residual=free_rms(trial_scaled,fixed)
              if(ieee_is_finite(trial_residual)) then
                if(.5_dp*trial_residual**2<=.5_dp*residual**2+options%armijo_c1*alpha*slope) then
                  accepted=.true.;exit
                end if
              end if
            end if
          end if
        end if
        alpha=alpha*options%backtrack_factor;report%total_backtracks=report%total_backtracks+1
      end do
      if(.not.accepted) exit
      update=alpha*free_rms(dx,fixed)
      state%x=trial%x;load=trial_load
    end do
    call set_status(status,STATUS_NONLINEAR_SOLVER_FAILED,'Temperature-control corrector failed; caller must roll back.')
  end subroutine solve_temperature_controlled_newton

  subroutine parameter_context(base,control,parameter,ctx,mapping,positions)
    type(SparseAssemblyContextType), intent(in) :: base
    integer, intent(in) :: control
    real(dp), intent(in) :: parameter(:)
    type(SparseAssemblyContextType), intent(out) :: ctx
    integer, allocatable, intent(out) :: mapping(:),positions(:)
    integer, allocatable :: columns(:)
    integer :: row,p,k
    logical :: inserted,needed
    ctx%number_of_equations=base%number_of_equations
    allocate(ctx%ia(ctx%number_of_equations+1),positions(ctx%number_of_equations), &
      mapping(base%number_of_nonzeros),columns(base%number_of_nonzeros+count(parameter/=0)))
    k=0;positions=0
    do row=1,ctx%number_of_equations
      ctx%ia(row)=k+1;inserted=.false.;needed=parameter(row)/=0
      do p=base%ia(row),base%ia(row+1)-1
        if(needed.and..not.inserted.and.base%ja(p)>control) then
          k=k+1;columns(k)=control;positions(row)=k;inserted=.true.
        end if
        k=k+1;columns(k)=base%ja(p);mapping(p)=k
        if(base%ja(p)==control) then
          positions(row)=k;inserted=.true.
        end if
      end do
      if(needed.and..not.inserted) then
        k=k+1;columns(k)=control;positions(row)=k
      end if
    end do
    ctx%ia(ctx%number_of_equations+1)=k+1;ctx%number_of_nonzeros=k
    allocate(ctx%ja(k));ctx%ja=columns(:k)
  end subroutine parameter_context

  real(dp) function free_rms(x,fixed) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: fixed(:)
    integer :: i,nfree
    value=0;nfree=0
    do i=1,size(x)
      if(fixed(i)) cycle
      value=value+x(i)*x(i);nfree=nfree+1
    end do
    value=sqrt(value/real(nfree,dp))
  end function free_rms
end module temperature_control_m

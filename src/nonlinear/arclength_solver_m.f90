! Opt-in pseudo-arclength diagnostic. Original PDE rows are retained.
! Lambda enters only the Dirichlet equations: their exact derivative is -db/dlambda.
! State coordinates and the continuation metric are scaled; the state metric is
! an RMS metric, so refining a mesh does not multiply the meaning of a step.
module arclength_solver_m
  use, intrinsic :: iso_fortran_env, only : output_unit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE,F_ETA_N,F_ETA_P,global_dof
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType,StateType,LoadStateType,NonlinearOptionsType, &
    ContinuationOptionsType,DiagnosticsType
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_INVALID_ARGUMENT,STATUS_CONTINUATION_FAILED
  use sparse_pattern_m, only : SparseAssemblyContextType
  use dirichlet_m, only : DirichletBCType,build_fixed_dof_mask,project_dirichlet_state
  use global_assembly_m, only : GlobalAssemblyRequestType,GlobalSystemType,assemble_global_system
  use scaling_m, only : build_variable_scaling,equilibrate_csr_system
  use pardiso_solver_m, only : PardisoDiagnosticsType,solve_csr_pardiso
  use newton_solver_m, only : NewtonReportType,solve_single_load_newton
  use balance_audit_m, only : audit_global_balances
  use continuation_solver_m, only : ContinuationReportType
  implicit none
  private
  public :: solve_arclength_continuation,trace_arclength_system

  abstract interface
    subroutine arc_evaluator(x,lambda,need_jacobian,residual,jacobian,parameter,status)
      import dp,SolverStatusType
      real(dp), intent(in) :: x(:),lambda
      logical, intent(in) :: need_jacobian
      real(dp), intent(out) :: residual(:),jacobian(:),parameter(:)
      type(SolverStatusType), intent(out) :: status
    end subroutine
    subroutine arc_checkpoint(state,load,point,residual,update,status)
      import dp,StateType,LoadStateType,SolverStatusType
      type(StateType), intent(in) :: state
      type(LoadStateType), intent(in) :: load
      integer, intent(in) :: point
      real(dp), intent(in) :: residual,update
      type(SolverStatusType), intent(out) :: status
    end subroutine
  end interface
contains
  subroutine solve_arclength_continuation(mesh,materials,context,bc,field_scales, &
      nonlinear_options,continuation_options,target,state,load,report,status,checkpoint)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: bc(:)
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    type(NonlinearOptionsType), intent(in) :: nonlinear_options
    type(ContinuationOptionsType), intent(in) :: continuation_options
    type(LoadStateType), intent(in) :: target
    type(StateType), intent(inout) :: state
    type(LoadStateType), intent(inout) :: load
    type(ContinuationReportType), intent(out) :: report
    type(SolverStatusType), intent(out) :: status
    procedure(arc_checkpoint), optional :: checkpoint
    real(dp), allocatable :: scales(:),parameter(:)
    logical, allocatable :: fixed(:),limited(:)
    type(NewtonReportType) :: nr
    type(StateType) :: eval_state,replay_state
    type(LoadStateType) :: eval_load,replay_load
    type(GlobalSystemType) :: system
    type(GlobalAssemblyRequestType) :: request
    type(DiagnosticsType) :: diag
    integer :: n,b,k,row
    call clear_status(status)
    report=ContinuationReportType()
    n=context%number_of_equations
    if(n/=NDOF_NODE*mesh%number_of_nodes().or.abs(target%joule_factor-load%joule_factor)>1e-12_dp) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Arclength requires matching mesh and a fixed Joule factor.');return
    end if
    if(nonlinear_options%trust_region) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Arclength uses bordered Armijo correction; disable newton_trust_region.')
      return
    end if
    call build_fixed_dof_mask(bc,mesh%number_of_nodes(),fixed,status)
    if(.not.status%is_ok()) return
    call build_variable_scaling(mesh%number_of_nodes(),field_scales,scales,status)
    if(.not.status%is_ok()) return
    allocate(parameter(n),limited(n),eval_state%x(n))
    parameter=0;limited=.false.
    limited(F_ETA_N::NDOF_NODE)=.true.;limited(F_ETA_P::NDOF_NODE)=.true.
    do b=1,size(bc)
      do k=1,size(bc(b)%node_ids)
        row=global_dof(bc(b)%node_ids(k),bc(b)%field_id)
        parameter(row)=-(bc(b)%value_end-bc(b)%value_start)
      end do
    end do
    call solve_single_load_newton(mesh,materials,load,context,bc,field_scales, &
      nonlinear_options,state,nr,status)
    if(.not.status%is_ok()) return
    eval_load=load
    call trace_arclength_system(context,scales,fixed,limited,nonlinear_options,continuation_options, &
      target%voltage_factor,evaluate,state,load,report,status,record)
    if(.not.status%is_ok()) return
    ! A crossed target is not an exact-voltage result: explicitly re-solve there.
    if(report%arclength_target_crossed.and..not.report%arclength_fold_detected) then
      replay_state=state;replay_load=load;replay_load%voltage_factor=target%voltage_factor
      call solve_single_load_newton(mesh,materials,replay_load,context,bc,field_scales, &
        nonlinear_options,replay_state,nr,status)
      if(.not.status%is_ok()) return
      call record(replay_state,replay_load,report%accepted_points+1,nr%final_residual_rms,nr%final_update_rms,status)
      if(.not.status%is_ok()) return
      k=report%accepted_points+1;report%accepted_points=k
      report%phase(k)=1;report%newton_iterations(k)=nr%iterations
      report%path_parameter(k)=report%path_parameter(k-1)
      report%voltage_factor(k)=replay_load%voltage_factor;report%joule_factor(k)=replay_load%joule_factor
      report%residual_rms(k)=nr%final_residual_rms;report%update_rms(k)=nr%final_update_rms
      report%total_newton_iterations=report%total_newton_iterations+nr%iterations
      state=replay_state;load=replay_load;report%converged=.true.;report%final_load=load
    end if
  contains
    subroutine evaluate(x,lambda,need_jacobian,residual,jacobian,dparameter,local_status)
      real(dp), intent(in) :: x(:),lambda
      logical, intent(in) :: need_jacobian
      real(dp), intent(out) :: residual(:),jacobian(:),dparameter(:)
      type(SolverStatusType), intent(out) :: local_status
      type(StateType) :: prescribed
      eval_state%x=x;eval_load%voltage_factor=lambda;request%need_jacobian=need_jacobian
      call assemble_global_system(mesh,materials,eval_state,eval_load,context,bc,request,system,diag,local_status)
      if(.not.local_status%is_ok()) return
      residual=system%residual;jacobian=system%jacobian_values;dparameter=parameter
      prescribed=eval_state
      call project_dirichlet_state(prescribed,bc,eval_load,mesh%number_of_nodes(),local_status)
      if(.not.local_status%is_ok()) return
      where(fixed) residual=x-prescribed%x
      ! Free rows contain J(x) including its boundary columns. Adding J*db/dlambda
      ! here would count the boundary lifting twice. Their partial dF/dlambda is zero.
    end subroutine
    subroutine record(saved,saved_load,index,residual,update,local_status)
      type(StateType), intent(in) :: saved
      type(LoadStateType), intent(in) :: saved_load
      integer, intent(in) :: index
      real(dp), intent(in) :: residual,update
      type(SolverStatusType), intent(out) :: local_status
      real(dp) :: balances(4),boundary_error,heat_free,heat_absolute
      type(DiagnosticsType) :: diagnostics
      call audit_global_balances(mesh,materials,saved,saved_load,context,bc,balances,boundary_error, &
        diagnostics,local_status,continuation_options%heat_balance_reference_power,heat_free,heat_absolute)
      if(.not.local_status%is_ok()) return
      if(maxval(balances)>=1e-8_dp.or.boundary_error>=1e-10_dp.or.diagnostics%min_temperature<=0.or. &
          diagnostics%min_n<=0.or.diagnostics%min_p<=0) then
        call set_status(local_status,STATUS_CONTINUATION_FAILED,'Arclength point failed algebraic/physical admissibility.');return
      end if
      ! The generic driver records residual/update; this callback supplies PDE audits.
      report%balance(:,index)=balances;report%boundary_error(index)=boundary_error
      report%heat_free_sum(index)=heat_free;report%heat_absolute_sum(index)=heat_absolute
      report%minimum_temperature(index)=diagnostics%min_temperature
      report%maximum_temperature(index)=diagnostics%max_temperature
      report%minimum_n(index)=diagnostics%min_n;report%maximum_n(index)=diagnostics%max_n
      report%minimum_p(index)=diagnostics%min_p;report%maximum_p(index)=diagnostics%max_p
      report%maximum_balance_defect=max(report%maximum_balance_defect,maxval(balances))
      report%maximum_boundary_error=max(report%maximum_boundary_error,boundary_error)
      if(present(checkpoint)) call checkpoint(saved,saved_load,index,residual,update,local_status)
    end subroutine
  end subroutine solve_arclength_continuation

  ! Generic backend also permits an analytic saddle-node regression using the
  ! same corrector, tangent solve, rollback and bracket refinement as the PDE.
  subroutine trace_arclength_system(base,scales,fixed,limited,options,arc_options,target,evaluate, &
      state,load,report,status,checkpoint)
    type(SparseAssemblyContextType), intent(in) :: base
    real(dp), intent(in) :: scales(:),target
    logical, intent(in) :: fixed(:),limited(:)
    type(NonlinearOptionsType), intent(in) :: options
    type(ContinuationOptionsType), intent(in) :: arc_options
    procedure(arc_evaluator) :: evaluate
    type(StateType), intent(inout) :: state
    type(LoadStateType), intent(inout) :: load
    type(ContinuationReportType), intent(out) :: report
    type(SolverStatusType), intent(out) :: status
    procedure(arc_checkpoint), optional :: checkpoint
    type(SparseAssemblyContextType) :: aug
    type(SolverStatusType) :: attempt_status
    type(StateType) :: saved
    type(LoadStateType) :: saved_load
    real(dp), allocatable :: z(:),trial(:),left(:),right(:),best(:),t(:),next_t(:),left_t(:),right_t(:)
    real(dp), allocatable :: weights(:),zscales(:),r(:),j(:),parameter(:),av(:),sv(:),sr(:),rs(:),dx(:),rhs(:)
    integer, allocatable :: mapping(:),border(:)
    real(dp) :: ds,residual,update,width,step,arclength,best_voltage
    integer :: n,i,p,pos,attempt,it,refinement,capacity
    call clear_status(status);report=ContinuationReportType()
    n=base%number_of_equations
    if(.not.arc_options%arclength_enabled.or.n<1.or.size(scales)/=n.or.size(state%x)/=n.or. &
        size(fixed)/=n.or.size(limited)/=n.or.count(.not.fixed)<1) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Invalid arclength dimensions/options.');return
    end if
    if(.not.all(ieee_is_finite([arc_options%arclength_psi_state,arc_options%arclength_psi_voltage, &
        arc_options%arclength_fold_tolerance,arc_options%initial_step,arc_options%minimum_step, &
        arc_options%maximum_step])).or.any(scales<=0).or..not.all(ieee_is_finite(scales)).or. &
        arc_options%arclength_psi_state<=0.or.arc_options%arclength_psi_voltage<=0.or. &
        arc_options%arclength_fold_tolerance<=0.or.arc_options%arclength_maximum_refinements<1.or. &
        arc_options%minimum_step<=0.or.arc_options%maximum_step<arc_options%initial_step.or. &
        arc_options%initial_step<=0.or.arc_options%maximum_attempts<1.or. &
        arc_options%maximum_newton_per_attempt<1.or.arc_options%cutback_factor<=0.or. &
        arc_options%cutback_factor>=1.or.options%minimum_step<=0.or.options%backtrack_factor<=0.or. &
        options%backtrack_factor>=1) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Arclength weights, step or refinement limits are invalid.');return
    end if
    allocate(z(n+1),trial(n+1),t(n+1),next_t(n+1),weights(n+1),zscales(n+1),r(n+1), &
      j(base%number_of_nonzeros),parameter(n),mapping(base%number_of_nonzeros),border(n),dx(n+1),rhs(n+1))
    weights(:n)=arc_options%arclength_psi_state/real(n,dp);weights(n+1)=arc_options%arclength_psi_voltage
    zscales(:n)=scales;zscales(n+1)=1
    z(:n)=state%x;z(n+1)=load%voltage_factor
    aug%number_of_equations=n+1;aug%number_of_nonzeros=base%number_of_nonzeros+2*n+1
    allocate(aug%ia(n+2),aug%ja(aug%number_of_nonzeros),av(aug%number_of_nonzeros))
    pos=1
    do i=1,n
      aug%ia(i)=pos
      do p=base%ia(i),base%ia(i+1)-1
        aug%ja(pos)=base%ja(p);mapping(p)=pos;pos=pos+1
      end do
      aug%ja(pos)=n+1;border(i)=pos;pos=pos+1
    end do
    aug%ia(n+1)=pos
    do i=1,n+1
      aug%ja(pos)=i;pos=pos+1
    end do
    aug%ia(n+2)=pos
    capacity=arc_options%maximum_attempts+arc_options%arclength_maximum_refinements+5
    call allocate_report(capacity)
    t=0;t(n+1)=1/sqrt(weights(n+1));arclength=0
    call assemble(z,z,t,0.0_dp,.true.,status)
    if(.not.status%is_ok()) return
    rhs=-r
    call equilibrate_csr_system(aug,av,rhs,zscales,sv,sr,rs,status)
    if(.not.status%is_ok()) return
    residual=free_rms(sr(:n))
    if(residual>options%residual_tolerance) then
      call set_status(status,STATUS_CONTINUATION_FAILED,'Arclength requires a converged initial state.');return
    end if
    if(any(fixed)) then
      if(maxval(abs(pack(r(:n)/scales,fixed)))>options%update_tolerance) then
        call set_status(status,STATUS_CONTINUATION_FAILED,'Arclength initial boundary mismatch.');return
      end if
    end if
    call save_point(z,0,residual,0.0_dp,status)
    if(.not.status%is_ok()) return
    if(z(n+1)>=target) then
      report%arclength_target_crossed=.true.;return
    end if
    call tangent(z,t,next_t,status)
    if(.not.status%is_ok()) return
    t=next_t;ds=arc_options%initial_step
    do attempt=1,arc_options%maximum_attempts
      report%total_attempts=attempt
      call correct(z,t,ds,trial,it,residual,update,attempt_status)
      if(.not.attempt_status%is_ok()) then
        report%rejected_attempts=report%rejected_attempts+1
        ds=ds*arc_options%cutback_factor
        if(options%verbose) then
          write(output_unit,'(a,es14.6,a,a)') 'arc_cutback=',ds,',reason=',trim(attempt_status%message)
          flush(output_unit)
        end if
        if(ds<arc_options%minimum_step) then
          status=attempt_status;return
        end if
        cycle
      end if
      call tangent(trial,t,next_t,status)
      if(.not.status%is_ok()) return
      arclength=arclength+ds
      call save_point(trial,it,residual,update,status)
      if(.not.status%is_ok()) return
      if(t(n+1)>0.and.next_t(n+1)<=0) then
        left=z;left_t=t;right=trial;right_t=next_t
        best=z;if(trial(n+1)>z(n+1)) best=trial
        best_voltage=best(n+1)
        report%arclength_fold_detected=.true.
        do refinement=1,arc_options%arclength_maximum_refinements
          width=arc_distance(left,right)
          if(width<=arc_options%arclength_fold_tolerance) exit
          step=0.5_dp*arc_projection(left_t,left,right)
          if(step<=0) then
            call set_status(status,STATUS_CONTINUATION_FAILED,'Fold bracket lost its orientation.');return
          end if
          call correct(left,left_t,step,trial,it,residual,update,status)
          if(.not.status%is_ok()) return
          call tangent(trial,left_t,next_t,status)
          if(.not.status%is_ok()) return
          call save_point(trial,it,residual,update,status)
          if(.not.status%is_ok()) return
          if(trial(n+1)>best_voltage) then
            best=trial;best_voltage=trial(n+1)
          end if
          if(next_t(n+1)>0) then
            left=trial;left_t=next_t
          else
            right=trial;right_t=next_t
          end if
        end do
        width=arc_distance(left,right)
        report%arclength_fold_bracket_width=width
        report%arclength_fold_left_voltage_factor=left(n+1)
        report%arclength_fold_right_voltage_factor=right(n+1)
        report%arclength_fold_left_tangent=left_t(n+1)
        report%arclength_fold_right_tangent=right_t(n+1)
        report%arclength_fold_voltage_factor=best_voltage
        report%arclength_fold_joule_factor=load%joule_factor
        report%arclength_fold_located=width<=arc_options%arclength_fold_tolerance
        state%x=best(:n);load%voltage_factor=best(n+1);report%final_load=load
        if(.not.report%arclength_fold_located) &
          call set_status(status,STATUS_CONTINUATION_FAILED,'Fold bracket did not reach requested arclength tolerance.')
        return
      end if
      z=trial;t=next_t
      state%x=z(:n);load%voltage_factor=z(n+1);report%final_load=load
      if(z(n+1)>=target) then
        report%arclength_target_crossed=.true.;return
      end if
      report%minimum_accepted_step=min(report%minimum_accepted_step,ds)
      if(it<=arc_options%fast_iteration_threshold) ds=min(arc_options%maximum_step,ds*arc_options%growth_factor)
      if(it>=arc_options%slow_iteration_threshold) ds=max(arc_options%minimum_step,ds*arc_options%cutback_factor)
    end do
    call set_status(status,STATUS_CONTINUATION_FAILED,'Arclength attempt limit; last accepted state retained.')
  contains
    subroutine assemble(point,anchor,direction,length,need_jacobian,local_status)
      real(dp), intent(in) :: point(:),anchor(:),direction(:),length
      logical, intent(in) :: need_jacobian
      type(SolverStatusType), intent(out) :: local_status
      call evaluate(point(:n),point(n+1),need_jacobian,r(:n),j,parameter,local_status)
      if(.not.local_status%is_ok()) return
      r(n+1)=0.0_dp
      do i=1,n+1
        r(n+1)=r(n+1)+weights(i)*direction(i)*(point(i)-anchor(i))/zscales(i)
      end do
      r(n+1)=r(n+1)-length
      if(need_jacobian) then
        av(mapping)=j;av(border)=parameter
        do pos=aug%ia(n+1),aug%ia(n+2)-1
          av(pos)=weights(pos-aug%ia(n+1)+1)*direction(pos-aug%ia(n+1)+1)/ &
            zscales(pos-aug%ia(n+1)+1)
        end do
      end if
      if(.not.all(ieee_is_finite(r))) &
        call set_status(local_status,STATUS_CONTINUATION_FAILED,'Non-finite arclength residual.')
    end subroutine
    subroutine linear_solve(rhs,solution,local_status)
      real(dp), intent(in) :: rhs(:)
      real(dp), intent(out) :: solution(:)
      type(SolverStatusType), intent(out) :: local_status
      type(PardisoDiagnosticsType) :: diagnostics
      call equilibrate_csr_system(aug,av,rhs,zscales,sv,sr,rs,local_status)
      if(.not.local_status%is_ok()) return
      call solve_csr_pardiso(n+1,aug%ia,aug%ja,sv,sr,solution,diagnostics,local_status,try_hole_split=.false.)
      ! Return dimensionless coordinates: tangent normalization and the update
      ! norm use y. Only predictor/candidate field updates multiply by zscales.
    end subroutine
    subroutine tangent(point,orientation,result,local_status)
      real(dp), intent(in) :: point(:),orientation(:)
      real(dp), intent(out) :: result(:)
      type(SolverStatusType), intent(out) :: local_status
      real(dp) :: norm
      call assemble(point,point,orientation,0.0_dp,.true.,local_status)
      if(.not.local_status%is_ok()) return
      rhs=0;rhs(n+1)=1
      call linear_solve(rhs,result,local_status)
      if(.not.local_status%is_ok()) return
      norm=arc_norm(result)
      if(.not.ieee_is_finite(norm).or.norm<=tiny(1.0_dp)) then
        call set_status(local_status,STATUS_CONTINUATION_FAILED,'Invalid branch tangent.');return
      end if
      do i=1,n+1
        result(i)=result(i)/norm
      end do
      if(arc_scaled_dot(result,orientation)<0) then
        do i=1,n+1
          result(i)=-result(i)
        end do
      end if
    end subroutine
    subroutine correct(anchor,direction,length,result,iterations,residual,update,local_status)
      real(dp), intent(in) :: anchor(:),direction(:),length
      real(dp), intent(out) :: result(:),residual,update
      integer, intent(out) :: iterations
      type(SolverStatusType), intent(out) :: local_status
      type(SolverStatusType) :: trial_status
      real(dp) :: base_merit,trial_merit,alpha,limit,logmax,boundary
      real(dp), allocatable :: candidate(:),frozen_scales(:)
      integer :: q
      logical :: accepted
      call clear_status(local_status)
      allocate(candidate(n+1),frozen_scales(n+1))
      do q=1,n+1
        result(q)=anchor(q)+length*direction(q)*zscales(q)
      end do
      update=huge(1.0_dp)
      do iterations=0,arc_options%maximum_newton_per_attempt
        call assemble(result,anchor,direction,length,.true.,local_status)
        if(.not.local_status%is_ok()) return
        rhs=-r
        call equilibrate_csr_system(aug,av,rhs,zscales,sv,sr,rs,local_status)
        if(.not.local_status%is_ok()) return
        residual=free_rms(sr(:n));boundary=0
        if(any(fixed)) then
          boundary=0.0_dp
          do q=1,n
            if(fixed(q)) boundary=max(boundary,abs(r(q)/scales(q)))
          end do
        end if
        if(options%verbose) then
          write(output_unit,'(a,i0,a,4(es14.6,a))') 'arc_iteration=',iterations,',lambda=',result(n+1), &
            ',residual=',residual,',update=',update,',constraint=',r(n+1),','
          flush(output_unit)
        end if
        if(residual<=options%residual_tolerance.and.abs(r(n+1))<=options%update_tolerance.and. &
            boundary<=options%update_tolerance.and.iterations>0.and.update<=options%update_tolerance) then
          return
        end if
        if(iterations==arc_options%maximum_newton_per_attempt) exit
        base_merit=0.5_dp*(free_rms(sr(:n))**2+r(n+1)**2)
        do q=1,n+1
          rhs(q)=-r(q)
        end do
        call linear_solve(rhs,dx,local_status)
        if(.not.local_status%is_ok()) return
        frozen_scales=rs;limit=1
        if(any(limited)) then
          logmax=0.0_dp
          do q=1,n
            if(limited(q)) logmax=max(logmax,abs(dx(q)*scales(q)))
          end do
          limit=min(1.0_dp,options%maximum_log_density_update/max(logmax,tiny(1.0_dp)))
        end if
        do q=1,n+1
          dx(q)=limit*dx(q)
        end do
        alpha=1;accepted=.false.
        do while(alpha>=options%minimum_step)
          do q=1,n+1
            candidate(q)=result(q)+alpha*dx(q)*zscales(q)
          end do
          if(all(ieee_is_finite(candidate)).and.candidate(n+1)>=0) then
            call assemble(candidate,anchor,direction,length,.false.,trial_status)
            if(trial_status%is_ok()) then
              trial_merit=0.5_dp*(scaled_free_rms(r(:n),frozen_scales(:n))**2+r(n+1)**2)
              if(trial_merit<=base_merit*(1-2*options%armijo_c1*alpha*limit)) then
                accepted=.true.;exit
              end if
            end if
          end if
          alpha=alpha*options%backtrack_factor
        end do
        if(.not.accepted) exit
        update=alpha*max(free_rms(dx(:n)),abs(dx(n+1)))
        result=candidate
      end do
      call set_status(local_status,STATUS_CONTINUATION_FAILED,'Bordered corrector failed; caller must roll back.')
    end subroutine
    real(dp) function free_rms(v) result(value)
      real(dp), intent(in) :: v(:)
      integer :: q
      value=0.0_dp
      do q=1,n
        if(.not.fixed(q)) value=value+v(q)*v(q)
      end do
      value=sqrt(value/real(count(.not.fixed),dp))
    end function
    real(dp) function scaled_free_rms(v,s) result(value)
      real(dp), intent(in) :: v(:),s(:)
      integer :: q
      value=0.0_dp
      do q=1,n
        if(.not.fixed(q)) value=value+(v(q)/s(q))**2
      end do
      value=sqrt(value/real(count(.not.fixed),dp))
    end function
    real(dp) function arc_dot(a,b) result(value)
      real(dp), intent(in) :: a(:),b(:)
      integer :: q
      value=0.0_dp
      do q=1,n+1
        value=value+weights(q)*a(q)*b(q)/zscales(q)
      end do
    end function
    real(dp) function arc_projection(a,left_point,right_point) result(value)
      real(dp), intent(in) :: a(:),left_point(:),right_point(:)
      integer :: q
      value=0.0_dp
      do q=1,n+1
        value=value+weights(q)*a(q)*(right_point(q)-left_point(q))/zscales(q)
      end do
    end function
    real(dp) function arc_scaled_dot(a,b) result(value)
      real(dp), intent(in) :: a(:),b(:)
      integer :: q
      value=0.0_dp
      do q=1,n+1
        value=value+weights(q)*a(q)*b(q)
      end do
    end function
    real(dp) function arc_norm(a) result(value)
      real(dp), intent(in) :: a(:)
      integer :: q
      value=0.0_dp
      do q=1,n+1
        value=value+weights(q)*a(q)*a(q)
      end do
      value=sqrt(value)
    end function
    real(dp) function arc_distance(a,b) result(value)
      real(dp), intent(in) :: a(:),b(:)
      integer :: q
      value=0.0_dp
      do q=1,n+1
        value=value+weights(q)*((a(q)-b(q))/zscales(q))**2
      end do
      value=sqrt(value)
    end function
    subroutine allocate_report(count)
      integer, intent(in) :: count
      allocate(report%phase(count),report%newton_iterations(count),report%path_parameter(count), &
        report%voltage_factor(count),report%joule_factor(count),report%residual_rms(count), &
        report%update_rms(count),report%balance(4,count),report%heat_free_sum(count),report%heat_absolute_sum(count), &
        report%boundary_error(count),report%minimum_temperature(count),report%maximum_temperature(count), &
        report%minimum_n(count),report%maximum_n(count),report%minimum_p(count),report%maximum_p(count))
      report%balance=0;report%boundary_error=0;report%heat_free_sum=0;report%heat_absolute_sum=0
    end subroutine
    subroutine save_point(point,iterations,residual,update,local_status)
      real(dp), intent(in) :: point(:),residual,update
      integer, intent(in) :: iterations
      type(SolverStatusType), intent(out) :: local_status
      integer :: k
      call clear_status(local_status)
      k=report%accepted_points+1;saved%x=point(:n);saved_load=load;saved_load%voltage_factor=point(n+1)
      if(present(checkpoint)) call checkpoint(saved,saved_load,k,residual,update,local_status)
      if(.not.local_status%is_ok()) return
      report%accepted_points=k;report%phase(k)=1
      report%voltage_factor(k)=point(n+1);report%joule_factor(k)=load%joule_factor
      report%newton_iterations(k)=iterations;report%residual_rms(k)=residual;report%update_rms(k)=update
      report%path_parameter(k)=arclength
      report%total_newton_iterations=report%total_newton_iterations+iterations
      if(options%verbose) then
        write(output_unit,'(a,i0,a,es24.16e3)') 'arc_accepted=',k,',voltage_factor=',point(n+1)
        flush(output_unit)
      end if
    end subroutine
  end subroutine trace_arclength_system
end module arclength_solver_m

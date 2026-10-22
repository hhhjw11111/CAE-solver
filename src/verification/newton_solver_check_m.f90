module newton_solver_check_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE, F_U1, F_U2, F_U3, F_PHI, F_ETA_N, &
    F_ETA_P, F_THETA, global_dof
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, &
    NonlinearOptionsType, DiagnosticsType
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED, &
    STATUS_NONLINEAR_SOLVER_FAILED
  use equilibrium_m, only : compute_equilibrium_densities
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context
  use dirichlet_m, only : DirichletBCType, project_dirichlet_state, build_fixed_dof_mask
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, &
    assemble_global_system
  use scaling_m, only : build_variable_scaling
  use newton_solver_m, only : NewtonReportType, solve_single_load_newton,dogleg_direction
  use material_tangent_check_m, only : define_test_material
  use global_assembly_check_m, only : build_two_element_mesh
  implicit none
  private

  real(dp), parameter :: FIELD_SCALES(NDOF_NODE)=[1.0e-10_dp,1.0e-10_dp, &
    1.0e-10_dp,10.0_dp,1.0e-1_dp,1.0e-1_dp,10.0_dp]

  public :: verify_newton_solver_g07

contains

  subroutine verify_newton_solver_g07(status,stagnation_check,trust_check)
    type(SolverStatusType), intent(out) :: status
    logical, intent(in), optional :: stagnation_check
    logical, intent(in), optional :: trust_check
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: warm_state,cold_state
    type(LoadStateType) :: load
    type(NonlinearOptionsType) :: options
    type(NewtonReportType) :: warm_report,cold_report
    type(SparseAssemblyContextType) :: context
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(SolverStatusType) :: local_status
    real(dp) :: n0,p0,solution_difference,boundary_error,warm_boundary_error,cold_boundary_error
    real(dp) :: warm_balance(4),cold_balance(4),balance(4)
    real(dp) :: minimum_temperature,maximum_temperature,minimum_n,maximum_n,minimum_p,maximum_p
    logical :: passed

    call clear_status(status)
    call build_two_element_mesh(mesh)
    mesh%coordinates=1.0e-6_dp*mesh%coordinates
    call define_test_material(materials(1))
    call compute_equilibrium_densities(materials(1)%donor_concentration, &
      materials(1)%acceptor_concentration,materials(1)%intrinsic_concentration, &
      n0,p0,local_status)
    if (.not. local_status%is_ok()) then
      call copy_verification_failure('G07 equilibrium state failed: ',local_status,status); return
    end if
    materials(1)%reference_n=n0
    materials(1)%reference_p=p0
    load%voltage_factor=1.0_dp
    load%joule_factor=1.0_dp
    options%max_newton=25
    options%residual_tolerance=1.0e-9_dp
    options%update_tolerance=1.0e-9_dp
    options%armijo_c1=1.0e-4_dp
    options%backtrack_factor=0.5_dp
    options%minimum_step=1.0e-8_dp
    options%maximum_log_density_update=1.0_dp
    if(present(stagnation_check)) then
      if(stagnation_check) then
        options%stagnation_window=4
      end if
    end if

    call define_g07_boundary_conditions(mesh,boundary_conditions,local_status)
    if (.not. local_status%is_ok()) then
      call copy_verification_failure('G07 boundary construction failed: ',local_status,status); return
    end if
    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call copy_verification_failure('G07 sparse initialization failed: ',local_status,status); return
    end if
    call define_initial_state(mesh,.true.,warm_state)
    call solve_single_load_newton(mesh,materials,load,context,boundary_conditions, &
      FIELD_SCALES,options,warm_state,warm_report,local_status)
    if (.not. local_status%is_ok()) then
      write(*,'(a)') 'verification_id=G07_true_residual_newton'
      call print_report('warm',warm_report)
      call copy_verification_failure('G07 warm-start Newton failed: ',local_status,status); return
    end if
    call define_initial_state(mesh,.false.,cold_state)
    if(present(trust_check)) then
      if(trust_check) then
        options%trust_region=.true.
        options%trust_initial_radius=0.025_dp
        options%max_newton=100
      end if
    end if
    call solve_single_load_newton(mesh,materials,load,context,boundary_conditions, &
      FIELD_SCALES,options,cold_state,cold_report,local_status)
    if (.not. local_status%is_ok()) then
      write(*,'(a)') 'verification_id=G07_true_residual_newton'
      call print_report('warm',warm_report)
      call print_report('cold',cold_report)
      call copy_verification_failure('G07 cold-start Newton failed: ',local_status,status); return
    end if

    solution_difference=scaled_state_difference(mesh,warm_state,cold_state,local_status)
    if (.not. local_status%is_ok()) then
      call copy_verification_failure('G07 solution comparison failed: ',local_status,status); return
    end if
    call audit_final_state(mesh,materials,load,context,boundary_conditions,warm_state, &
      warm_balance,warm_boundary_error,local_status)
    if (.not. local_status%is_ok()) then
      call copy_verification_failure('G07 warm-state audit failed: ',local_status,status); return
    end if
    call audit_final_state(mesh,materials,load,context,boundary_conditions,cold_state, &
      cold_balance,cold_boundary_error,local_status)
    if (.not. local_status%is_ok()) then
      call copy_verification_failure('G07 cold-state audit failed: ',local_status,status); return
    end if
    balance=max(warm_balance,cold_balance)
    boundary_error=max(warm_boundary_error,cold_boundary_error)
    minimum_temperature=min(warm_report%final_diagnostics%min_temperature, &
      cold_report%final_diagnostics%min_temperature)
    maximum_temperature=max(warm_report%final_diagnostics%max_temperature, &
      cold_report%final_diagnostics%max_temperature)
    minimum_n=min(warm_report%final_diagnostics%min_n,cold_report%final_diagnostics%min_n)
    maximum_n=max(warm_report%final_diagnostics%max_n,cold_report%final_diagnostics%max_n)
    minimum_p=min(warm_report%final_diagnostics%min_p,cold_report%final_diagnostics%min_p)
    maximum_p=max(warm_report%final_diagnostics%max_p,cold_report%final_diagnostics%max_p)

    passed=warm_report%converged .and. cold_report%converged .and. &
      warm_report%iterations <= options%max_newton .and. &
      cold_report%iterations <= options%max_newton .and. &
      warm_report%final_residual_rms < options%residual_tolerance .and. &
      cold_report%final_residual_rms < options%residual_tolerance .and. &
      warm_report%final_update_rms < options%update_tolerance .and. &
      cold_report%final_update_rms < options%update_tolerance .and. &
      warm_report%maximum_armijo_violation <= 1.0e-14_dp .and. &
      cold_report%maximum_armijo_violation <= 1.0e-14_dp .and. &
      solution_difference < 1.0e-7_dp .and. maxval(balance) < 1.0e-8_dp .and. &
      boundary_error < 1.0e-14_dp .and. minimum_temperature > 0.0_dp .and. &
      minimum_n > 0.0_dp .and. minimum_p > 0.0_dp .and. &
      ieee_is_finite(minimum_temperature) .and. ieee_is_finite(maximum_temperature) .and. &
      ieee_is_finite(minimum_n) .and. ieee_is_finite(maximum_n) .and. &
      ieee_is_finite(minimum_p) .and. ieee_is_finite(maximum_p)

    if(present(stagnation_check)) then
      if(stagnation_check) then
        block
          type(StateType) :: stalled_state
          type(NewtonReportType) :: stalled_report
          call define_initial_state(mesh,.false.,stalled_state)
          options%global_log_step_limit=.true.
          options%maximum_log_density_update=1.0e-8_dp
          call solve_single_load_newton(mesh,materials,load,context,boundary_conditions, &
            FIELD_SCALES,options,stalled_state,stalled_report,local_status)
          passed=passed .and. local_status%code==STATUS_NONLINEAR_SOLVER_FAILED .and. &
            index(local_status%message,'stagnation')>0 .and. .not.stalled_report%converged .and. &
            stalled_report%iterations==options%stagnation_window
          write(*,'(a,i0)') 'g07_stagnation_rejected_after_iterations=',stalled_report%iterations
          write(*,'(a)') 'g07_stagnation_message='//trim(local_status%message)
        end block
      end if
    end if

    if(present(trust_check)) then
      if(trust_check) then
        block
          real(dp) :: pn(3),pc(3),direction(3),radius,norm_d
          logical :: mask(3)
          integer :: test
          type(StateType) :: rejected_state
          type(NewtonReportType) :: rejected_report
          pn=[4.0_dp,2.0_dp,99.0_dp];pc=[0.5_dp,0.5_dp,0.0_dp];mask=[.false.,.false.,.true.]
          do test=1,3
            select case(test)
            case(1);radius=0.1_dp
            case(2);radius=1.0_dp
            case(3);radius=5.0_dp
            end select
            call dogleg_direction(pn,pc,mask,radius,direction)
            norm_d=sqrt(sum(direction(1:2)**2)/2.0_dp)
            passed=passed.and.direction(3)==0.0_dp.and.norm_d<=radius*(1.0_dp+1.0e-14_dp)
            if(test<3) passed=passed.and.abs(norm_d-radius)<1.0e-14_dp
            if(test==3) passed=passed.and.maxval(abs(direction(1:2)-pn(1:2)))==0.0_dp
          end do
          passed=passed.and.cold_report%trust_region_steps>0.and.warm_report%trust_region_steps==0
          call define_initial_state(mesh,.false.,rejected_state)
          options%trust_initial_radius=-1.0_dp
          call solve_single_load_newton(mesh,materials,load,context,boundary_conditions, &
            FIELD_SCALES,options,rejected_state,rejected_report,local_status)
          passed=passed.and..not.local_status%is_ok().and..not.rejected_report%converged
          write(*,'(a,i0)') 'g07_trust_accepted_steps=',cold_report%trust_region_steps
          write(*,'(a,i0)') 'g07_trust_rejected_trials=',cold_report%trust_region_rejections
          write(*,'(a)') 'g07_trust_comparison=LEGACY_WARM_VS_TRUST_COLD'
        end block
      end if
    end if

    write(*,'(a)') 'verification_id=G07_true_residual_newton'
    call print_report('warm',warm_report)
    call print_report('cold',cold_report)
    call print_metric('g07_scaled_solution_difference',solution_difference)
    call print_metric('g07_charge_balance_defect',balance(1))
    call print_metric('g07_electron_current_balance_defect',balance(2))
    call print_metric('g07_hole_current_balance_defect',balance(3))
    call print_metric('g07_heat_balance_defect',balance(4))
    call print_metric('g07_boundary_max_error',boundary_error)
    call print_metric('g07_minimum_temperature',minimum_temperature)
    call print_metric('g07_maximum_temperature',maximum_temperature)
    call print_metric('g07_minimum_electron_concentration',minimum_n)
    call print_metric('g07_maximum_electron_concentration',maximum_n)
    call print_metric('g07_minimum_hole_concentration',minimum_p)
    call print_metric('g07_maximum_hole_concentration',maximum_p)
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'G07 Newton convergence, Armijo, agreement, conservation, or admissibility failed.')
    end if
  end subroutine verify_newton_solver_g07

  subroutine define_g07_boundary_conditions(mesh,boundary_conditions,status)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), allocatable, intent(out) :: boundary_conditions(:)
    type(SolverStatusType), intent(out) :: status
    integer, allocatable :: left_nodes(:),right_nodes(:)

    call clear_status(status)
    call collect_face_nodes(mesh,0.0_dp,left_nodes)
    call collect_face_nodes(mesh,2.0e-6_dp,right_nodes)
    if (size(left_nodes) == 0 .or. size(right_nodes) == 0) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'G07 contact face is empty.')
      return
    end if
    allocate(boundary_conditions(11))
    call set_boundary(boundary_conditions(1),F_U1,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(2),F_U2,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(3),F_U3,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(4),F_PHI,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(5),F_PHI,right_nodes,0.0_dp,-10.0_dp)
    call set_boundary(boundary_conditions(6),F_ETA_N,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(7),F_ETA_N,right_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(8),F_ETA_P,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(9),F_ETA_P,right_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(10),F_THETA,left_nodes,0.0_dp,0.0_dp)
    call set_boundary(boundary_conditions(11),F_THETA,right_nodes,0.0_dp,0.0_dp)
  end subroutine define_g07_boundary_conditions

  subroutine collect_face_nodes(mesh,x_coordinate,nodes)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: x_coordinate
    integer, allocatable, intent(out) :: nodes(:)
    integer :: node,entry

    allocate(nodes(count(abs(mesh%coordinates(1,:)-x_coordinate) <= 1.0e-16_dp)))
    entry=0
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(1,node)-x_coordinate) > 1.0e-16_dp) cycle
      entry=entry+1
      nodes(entry)=node
    end do
  end subroutine collect_face_nodes

  subroutine set_boundary(boundary,field,nodes,value_start,value_end)
    type(DirichletBCType), intent(out) :: boundary
    integer, intent(in) :: field,nodes(:)
    real(dp), intent(in) :: value_start,value_end

    boundary%field_id=field
    allocate(boundary%node_ids(size(nodes)))
    boundary%node_ids=nodes
    boundary%value_start=value_start
    boundary%value_end=value_end
  end subroutine set_boundary

  subroutine define_initial_state(mesh,warm_start,state)
    type(MeshType), intent(in) :: mesh
    logical, intent(in) :: warm_start
    type(StateType), intent(out) :: state
    integer :: node

    allocate(state%x(NDOF_NODE*mesh%number_of_nodes()))
    state%x=0.0_dp
    if (.not. warm_start) return
    do node=1,mesh%number_of_nodes()
      state%x(global_dof(node,F_PHI))=-10.0_dp*mesh%coordinates(1,node)/2.0e-6_dp
    end do
  end subroutine define_initial_state

  real(dp) function scaled_state_difference(mesh,state_a,state_b,status) result(value)
    type(MeshType), intent(in) :: mesh
    type(StateType), intent(in) :: state_a,state_b
    type(SolverStatusType), intent(out) :: status
    real(dp), allocatable :: variable_scales(:)
    integer :: equation

    call build_variable_scaling(mesh%number_of_nodes(),FIELD_SCALES,variable_scales,status)
    if (.not. status%is_ok()) then
      value=huge(1.0_dp); return
    end if
    value=0.0_dp
    do equation=1,size(state_a%x)
      value=value+((state_a%x(equation)-state_b%x(equation))/variable_scales(equation))**2
    end do
    value=sqrt(value/real(size(state_a%x),dp))
  end function scaled_state_difference

  subroutine audit_final_state(mesh,materials,load,context,boundary_conditions,state, &
                               balance,boundary_error,status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(LoadStateType), intent(in) :: load
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    type(StateType), intent(in) :: state
    real(dp), intent(out) :: balance(4)
    real(dp), intent(out) :: boundary_error
    type(SolverStatusType), intent(out) :: status
    type(StateType) :: projected_state
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: unconstrained_system
    type(DiagnosticsType) :: diagnostics
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
  end subroutine audit_final_state

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

  subroutine print_report(label,report)
    character(len=*), intent(in) :: label
    type(NewtonReportType), intent(in) :: report
    integer :: entry

    write(*,'(a,l1)') 'g07_'//trim(label)//'_converged=',report%converged
    write(*,'(a,i0)') 'g07_'//trim(label)//'_iterations=',report%iterations
    write(*,'(a,i0)') 'g07_'//trim(label)//'_linear_solves=',report%total_linear_solves
    write(*,'(a,i0)') 'g07_'//trim(label)//'_backtracks=',report%total_backtracks
    write(*,'(a,i0)') 'g07_'//trim(label)//'_log_limited_steps=',report%log_limited_steps
    write(*,'(a,i0)') 'g07_'//trim(label)//'_maximum_perturbed_pivots=', &
      report%maximum_perturbed_pivots
    call print_metric('g07_'//trim(label)//'_initial_residual_rms', &
      report%initial_residual_rms)
    call print_metric('g07_'//trim(label)//'_final_residual_rms', &
      report%final_residual_rms)
    call print_metric('g07_'//trim(label)//'_final_update_rms',report%final_update_rms)
    call print_metric('g07_'//trim(label)//'_minimum_accepted_step', &
      report%minimum_accepted_step)
    call print_metric('g07_'//trim(label)//'_maximum_armijo_violation', &
      report%maximum_armijo_violation)
    call print_metric('g07_'//trim(label)//'_final_update_u1',report%final_update_by_field(F_U1))
    call print_metric('g07_'//trim(label)//'_final_update_u2',report%final_update_by_field(F_U2))
    call print_metric('g07_'//trim(label)//'_final_update_u3',report%final_update_by_field(F_U3))
    call print_metric('g07_'//trim(label)//'_final_update_phi',report%final_update_by_field(F_PHI))
    call print_metric('g07_'//trim(label)//'_final_update_eta_n', &
      report%final_update_by_field(F_ETA_N))
    call print_metric('g07_'//trim(label)//'_final_update_eta_p', &
      report%final_update_by_field(F_ETA_P))
    call print_metric('g07_'//trim(label)//'_final_update_theta', &
      report%final_update_by_field(F_THETA))
    call print_metric('g07_'//trim(label)//'_minimum_temperature', &
      report%final_diagnostics%min_temperature)
    call print_metric('g07_'//trim(label)//'_maximum_temperature', &
      report%final_diagnostics%max_temperature)
    call print_metric('g07_'//trim(label)//'_minimum_n',report%final_diagnostics%min_n)
    call print_metric('g07_'//trim(label)//'_maximum_n',report%final_diagnostics%max_n)
    call print_metric('g07_'//trim(label)//'_minimum_p',report%final_diagnostics%min_p)
    call print_metric('g07_'//trim(label)//'_maximum_p',report%final_diagnostics%max_p)
    do entry=1,report%history_count
      call print_indexed_metric('g07_'//trim(label)//'_residual',entry, &
        report%residual_history(entry))
    end do
    do entry=1,report%iterations
      call print_indexed_metric('g07_'//trim(label)//'_update',entry, &
        report%update_history(entry))
      call print_indexed_metric('g07_'//trim(label)//'_step',entry,report%step_history(entry))
    end do
  end subroutine print_report

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

  subroutine print_indexed_metric(name,index,value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: index
    real(dp), intent(in) :: value
    character(len=16) :: index_text

    write(index_text,'(i0)') index
    call print_metric(trim(name)//'_'//trim(index_text),value)
  end subroutine print_indexed_metric

  subroutine copy_verification_failure(prefix,source,destination)
    character(len=*), intent(in) :: prefix
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,STATUS_VERIFICATION_FAILED,trim(prefix)//trim(source%message))
  end subroutine copy_verification_failure

end module newton_solver_check_m

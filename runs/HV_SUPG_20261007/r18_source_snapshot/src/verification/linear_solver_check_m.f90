module linear_solver_check_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : NDOF_NODE, F_ETA_N, F_ETA_P, F_THETA
  use mesh_types_m, only : MeshType
  use model_types_m, only : MaterialType, StateType, LoadStateType, DiagnosticsType
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context, csr_matvec
  use dirichlet_m, only : DirichletBCType, project_dirichlet_state, build_fixed_dof_mask
  use global_assembly_m, only : GlobalAssemblyRequestType, GlobalSystemType, assemble_global_system
  use pardiso_solver_m, only : PardisoDiagnosticsType, solve_csr_pardiso
  use scaling_m, only : build_variable_scaling, equilibrate_csr_system, unscale_solution
  use material_tangent_check_m, only : define_test_material
  use global_assembly_check_m, only : build_two_element_mesh, define_global_state
  implicit none
  private

  public :: verify_linear_solver_l01,verify_hole_split_linear_solver

contains

  subroutine verify_hole_split_linear_solver(status)
    type(SolverStatusType), intent(out) :: status
    integer, parameter :: n=14,nnz=n*n
    integer :: ia(n+1),ja(nnz),row,col,k,example
    real(dp) :: a(n,n),exact(n),rhs(n),values(nnz),solution(n),full(n),error,full_error
    type(PardisoDiagnosticsType) :: diagnostics,full_diagnostics
    type(SolverStatusType) :: local_status
    logical :: passed,row_hole,col_hole
    call clear_status(status);passed=.true.
    write(*,'(a)') 'verification_id=BLK_checked_hole_block_linear_solve'
    do example=1,3
      a=0
      do row=1,n
        exact(row)=sin(.3_dp*row)+.1_dp
        row_hole=mod(row-1,NDOF_NODE)+1==F_ETA_P
        if(example==3.and.row_hole) exact(row)=1.0e14_dp
        do col=1,n
          col_hole=mod(col-1,NDOF_NODE)+1==F_ETA_P
          if(row==col) then
            a(row,col)=4.0_dp+.1_dp*row
          else if(row_hole.and..not.col_hole) then
            a(row,col)=.03_dp*sin(real(row+col,dp))
          else if(.not.row_hole.and.col_hole) then
            if(example==2) a(row,col)=.1_dp
            if(example==3) a(row,col)=1.0e-14_dp
          else if(abs(row-col)==1) then
            a(row,col)=-.2_dp
          end if
        end do
      end do
      ! Include a zero incremental Dirichlet row in the manufactured systems.
      a(1,:)=0;a(1,1)=1;exact(1)=0
      rhs=matmul(a,exact);k=0
      do row=1,n
        ia(row)=k+1
        do col=1,n
          k=k+1;ja(k)=col;values(k)=a(row,col)
        end do
      end do
      ia(n+1)=k+1
      call solve_csr_pardiso(n,ia,ja,values,rhs,full,full_diagnostics,local_status,try_hole_split=.false.)
      if(.not.local_status%is_ok()) then
        status=local_status;return
      end if
      call solve_csr_pardiso(n,ia,ja,values,rhs,solution,diagnostics,local_status,try_hole_split=.true.)
      if(.not.local_status%is_ok()) then
        status=local_status;return
      end if
      error=maxval(abs(solution-exact)/max(1.0_dp,abs(exact)))
      full_error=maxval(abs(solution-full)/max(1.0_dp,abs(full)))
      passed=passed.and.error<1.0e-12_dp.and.full_error<1.0e-12_dp.and.diagnostics%split_attempted
      if(example==1) passed=passed.and.diagnostics%split_accepted
      if(example>=2) passed=passed.and..not.diagnostics%split_accepted
      if(example==3) passed=passed.and.diagnostics%split_backward_error>1.0e-11_dp &
        .and.diagnostics%split_backward_error<1.0_dp
      write(*,'(a,i0,a,l1,a,es12.4,a,es12.4)') 'block_case=',example, &
        ',split_accepted=',diagnostics%split_accepted,',solution_error=',error, &
        ',full_matrix_difference=',full_error
    end do
    if(passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'Checked hole-block solve or full-matrix fallback failed.')
    end if
  end subroutine verify_hole_split_linear_solver

  subroutine verify_linear_solver_l01(status)
    type(SolverStatusType), intent(out) :: status
    type(SolverStatusType) :: local_status
    type(PardisoDiagnosticsType) :: small_diagnostics,global_diagnostics
    real(dp) :: small_solution_error,small_residual_error
    real(dp) :: global_solution_error,global_residual_error,scaled_matrix_max
    integer :: global_equations,global_nonzeros
    logical :: passed

    call clear_status(status)
    call verify_small_nonsymmetric_system(small_solution_error,small_residual_error, &
      small_diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call verify_scaled_global_system(global_solution_error,global_residual_error, &
      scaled_matrix_max,global_equations,global_nonzeros,global_diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if

    passed = small_solution_error < 1.0e-13_dp .and. small_residual_error < 1.0e-13_dp .and. &
      global_solution_error < 2.0e-8_dp .and. global_residual_error < 2.0e-10_dp .and. &
      scaled_matrix_max <= 1.0_dp+1.0e-14_dp .and. small_diagnostics%error_code == 0 .and. &
      global_diagnostics%error_code == 0

    write(*,'(a)') 'verification_id=L01_linear_solver'
    call print_metric('small_solution_relative_error',small_solution_error)
    call print_metric('small_residual_relative_error',small_residual_error)
    write(*,'(a,i0)') 'global_equation_count=',global_equations
    write(*,'(a,i0)') 'global_csr_nnz=',global_nonzeros
    call print_metric('scaled_global_matrix_max_abs',scaled_matrix_max)
    call print_metric('global_scaled_solution_relative_error',global_solution_error)
    call print_metric('global_physical_residual_relative_error',global_residual_error)
    write(*,'(a,i0)') 'pardiso_factor_nonzeros=',global_diagnostics%factor_nonzeros
    write(*,'(a,i0)') 'pardiso_perturbed_pivots=',global_diagnostics%perturbed_pivots
    write(*,'(a,i0)') 'pardiso_iterative_refinements=',global_diagnostics%iterative_refinements
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        'L01 PARDISO, scaling, or manufactured global solve exceeded tolerance.')
    end if
  end subroutine verify_linear_solver_l01

  subroutine verify_small_nonsymmetric_system(solution_error,residual_error,diagnostics,status)
    real(dp), intent(out) :: solution_error,residual_error
    type(PardisoDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    integer, parameter :: n=5,nnz=16
    integer :: ia(n+1),ja(nnz),row,position
    real(dp) :: values(nnz),exact(n),rhs(n),solution(n),product(n)

    ia=[1,4,7,11,14,17]
    ja=[1,2,5, 1,2,3, 2,3,4,5, 1,4,5, 2,4,5]
    values=[10.0_dp,2.0_dp,-1.0_dp, 3.0_dp,15.0_dp,4.0_dp, &
      -2.0_dp,12.0_dp,1.0_dp,1.0_dp, 1.0_dp,9.0_dp,2.0_dp, &
      2.0_dp,-3.0_dp,11.0_dp]
    exact=[1.0_dp,-2.0_dp,0.5_dp,3.0_dp,-1.0_dp]
    rhs=0.0_dp
    do row=1,n
      do position=ia(row),ia(row+1)-1
        rhs(row)=rhs(row)+values(position)*exact(ja(position))
      end do
    end do
    call solve_csr_pardiso(n,ia,ja,values,rhs,solution,diagnostics,status)
    if (.not. status%is_ok()) return
    product=0.0_dp
    do row=1,n
      do position=ia(row),ia(row+1)-1
        product(row)=product(row)+values(position)*solution(ja(position))
      end do
    end do
    solution_error=relative_norm(solution-exact,exact)
    residual_error=relative_norm(product-rhs,rhs)
  end subroutine verify_small_nonsymmetric_system

  subroutine verify_scaled_global_system(solution_error,residual_error,scaled_matrix_max, &
                                         equation_count,nonzero_count,diagnostics,status)
    real(dp), intent(out) :: solution_error,residual_error,scaled_matrix_max
    integer, intent(out) :: equation_count,nonzero_count
    type(PardisoDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(MaterialType) :: materials(1)
    type(StateType) :: state
    type(LoadStateType) :: load
    type(DiagnosticsType) :: assembly_diagnostics
    type(SparseAssemblyContextType) :: context
    type(GlobalAssemblyRequestType) :: request
    type(GlobalSystemType) :: system
    type(DirichletBCType), allocatable :: boundary_conditions(:)
    type(SolverStatusType) :: local_status
    logical, allocatable :: fixed(:)
    real(dp), allocatable :: exact(:),exact_scaled(:),rhs(:),solution(:),solution_scaled(:)
    real(dp), allocatable :: product(:),variable_scales(:),residual_scales(:)
    real(dp), allocatable :: scaled_values(:),scaled_rhs(:)
    real(dp) :: field_scales(NDOF_NODE)
    integer :: equation,node_count

    call clear_status(status)
    call build_two_element_mesh(mesh)
    call define_test_material(materials(1))
    materials(1)%reference_p=2.0e20_dp
    node_count=mesh%number_of_nodes()
    allocate(state%x(NDOF_NODE*node_count))
    call define_global_state(mesh,state)
    load%voltage_factor=0.60_dp
    load%joule_factor=0.77_dp
    call define_left_face_boundary_conditions(mesh,boundary_conditions)
    call project_dirichlet_state(state,boundary_conditions,load,node_count,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    request%need_jacobian=.true.
    call assemble_global_system(mesh,materials,state,load,context,boundary_conditions, &
      request,system,assembly_diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    call build_fixed_dof_mask(boundary_conditions,node_count,fixed,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if

    field_scales=[1.0e-4_dp,1.0e-4_dp,1.0e-4_dp,2.0e5_dp,1.0_dp,1.0_dp,1.0e2_dp]
    call build_variable_scaling(node_count,field_scales,variable_scales,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    equation_count=context%number_of_equations
    nonzero_count=context%number_of_nonzeros
    allocate(exact(equation_count),exact_scaled(equation_count),rhs(equation_count), &
      solution(equation_count),solution_scaled(equation_count),product(equation_count))
    do equation=1,equation_count
      exact_scaled(equation)=0.2_dp*sin(0.13_dp*real(equation,dp)) + &
        0.1_dp*cos(0.29_dp*real(equation,dp))
      if (fixed(equation)) exact_scaled(equation)=0.0_dp
    end do
    exact=variable_scales*exact_scaled
    call csr_matvec(context,system%jacobian_values,exact,rhs,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    call equilibrate_csr_system(context,system%jacobian_values,rhs,variable_scales, &
      scaled_values,scaled_rhs,residual_scales,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    scaled_matrix_max=maxval(abs(scaled_values))
    call solve_csr_pardiso(equation_count,context%ia,context%ja,scaled_values, &
      scaled_rhs,solution_scaled,diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    call unscale_solution(variable_scales,solution_scaled,solution,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    call csr_matvec(context,system%jacobian_values,solution,product,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message)); return
    end if
    solution_error=relative_norm(solution_scaled-exact_scaled,exact_scaled)
    residual_error=relative_norm(product-rhs,rhs)
  end subroutine verify_scaled_global_system

  subroutine define_left_face_boundary_conditions(mesh,boundary_conditions)
    type(MeshType), intent(in) :: mesh
    type(DirichletBCType), allocatable, intent(out) :: boundary_conditions(:)
    integer, allocatable :: left_nodes(:)
    integer :: node,face_count,field

    face_count=count(abs(mesh%coordinates(1,:)) < 1.0e-14_dp)
    allocate(left_nodes(face_count))
    face_count=0
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(1,node)) < 1.0e-14_dp) then
        face_count=face_count+1
        left_nodes(face_count)=node
      end if
    end do
    allocate(boundary_conditions(NDOF_NODE))
    do field=1,NDOF_NODE
      boundary_conditions(field)%field_id=field
      allocate(boundary_conditions(field)%node_ids(size(left_nodes)))
      boundary_conditions(field)%node_ids=left_nodes
    end do
    boundary_conditions(F_ETA_N)%value_start=0.12_dp
    boundary_conditions(F_ETA_N)%value_end=0.12_dp
    boundary_conditions(F_ETA_P)%value_start=-0.18_dp
    boundary_conditions(F_ETA_P)%value_end=-0.18_dp
    boundary_conditions(F_THETA)%value_start=45.0_dp
    boundary_conditions(F_THETA)%value_end=45.0_dp
  end subroutine define_left_face_boundary_conditions

  real(dp) function relative_norm(difference,reference) result(value)
    real(dp), intent(in) :: difference(:),reference(:)
    value=sqrt(sum(difference*difference))/max(sqrt(sum(reference*reference)),tiny(1.0_dp))
  end function relative_norm

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

end module linear_solver_check_m

module analytic_benchmarks_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : F_U1, F_PHI, F_ETA_N, F_ETA_P, F_THETA, global_dof
  use mesh_types_m, only : MeshType
  use pardiso_solver_m, only : PardisoDiagnosticsType
  use m0_benchmark_m, only : M0CoefficientsType, M0_V02, M0_V03, M0_V04, &
    M0_V05_LINEAR, M0_V05_UNIFORM, M0_V06, solve_m0_benchmark, &
    compute_m0_field_error, compute_m0_max_current, compute_m0_uniform_heat_balance
  implicit none
  private

  integer, parameter :: CONVERGENCE_LEVELS=5
  integer, parameter :: ELEMENT_COUNTS(CONVERGENCE_LEVELS)=[1,2,4,8,16]

  public :: run_m0_verification

contains

  subroutine run_m0_verification(verification_id,status)
    character(len=*), intent(in) :: verification_id
    type(SolverStatusType), intent(out) :: status
    type(SolverStatusType) :: local_status
    logical :: passed,current_pass

    call clear_status(status)
    passed=.true.
    write(*,'(a)') 'verification_id='//trim(uppercase(verification_id))//'_M0_analytic'
    select case(uppercase(trim(verification_id)))
    case('V02')
      call check_v02(current_pass,local_status)
      passed=current_pass
    case('V03')
      call check_v03(current_pass,local_status)
      passed=current_pass
    case('V04')
      call check_v04(current_pass,local_status)
      passed=current_pass
    case('V05')
      call check_v05(current_pass,local_status)
      passed=current_pass
    case('V06')
      call check_v06(current_pass,local_status)
      passed=current_pass
    case('G05')
      call check_v02(current_pass,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      passed=passed .and. current_pass
      call check_v03(current_pass,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      passed=passed .and. current_pass
      call check_v04(current_pass,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      passed=passed .and. current_pass
      call check_v05(current_pass,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      passed=passed .and. current_pass
      call check_v06(current_pass,local_status)
      if (.not. local_status%is_ok()) then
        call copy_failure(local_status,status); return
      end if
      passed=passed .and. current_pass
    case default
      call set_status(status,STATUS_VERIFICATION_FAILED,'Unknown M0 verification identifier.')
      return
    end select
    if (.not. local_status%is_ok()) then
      call copy_failure(local_status,status)
      return
    end if
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED, &
        trim(uppercase(verification_id))//' M0 analytic metric exceeded its frozen threshold.')
    end if
  end subroutine run_m0_verification

  subroutine check_v02(passed,status)
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(M0CoefficientsType) :: coefficients
    type(MeshType) :: mesh
    type(PardisoDiagnosticsType) :: diagnostics
    real(dp), allocatable :: solution(:)
    real(dp) :: phi_l2(CONVERGENCE_LEVELS),phi_h1(CONVERGENCE_LEVELS)
    real(dp) :: u_l2(CONVERGENCE_LEVELS),u_h1(CONVERGENCE_LEVELS)
    real(dp) :: phi_l2_order,phi_h1_order,u_l2_order,u_h1_order
    integer :: level,max_perturbed_pivots

    call clear_status(status)
    max_perturbed_pivots=0
    do level=1,CONVERGENCE_LEVELS
      call solve_m0_benchmark(ELEMENT_COUNTS(level),M0_V02,coefficients,mesh,solution, &
        diagnostics,status)
      if (.not. status%is_ok()) return
      max_perturbed_pivots=max(max_perturbed_pivots,diagnostics%perturbed_pivots)
      call compute_m0_field_error(mesh,solution,M0_V02,F_PHI,coefficients, &
        phi_l2(level),phi_h1(level),status)
      if (.not. status%is_ok()) return
      call compute_m0_field_error(mesh,solution,M0_V02,F_U1,coefficients, &
        u_l2(level),u_h1(level),status)
      if (.not. status%is_ok()) return
      call print_level_metrics('v02',ELEMENT_COUNTS(level),phi_l2(level),phi_h1(level))
    end do
    phi_l2_order=observed_order(phi_l2)
    phi_h1_order=observed_order(phi_h1)
    u_l2_order=observed_order(u_l2)
    u_h1_order=observed_order(u_h1)
    call print_metric('v02_phi_l2_order',phi_l2_order)
    call print_metric('v02_phi_h1_order',phi_h1_order)
    call print_metric('v02_u_l2_final',u_l2(CONVERGENCE_LEVELS))
    call print_metric('v02_u_h1_final',u_h1(CONVERGENCE_LEVELS))
    call print_metric('v02_u_l2_order',u_l2_order)
    call print_metric('v02_u_h1_order',u_h1_order)
    write(*,'(a,i0)') 'v02_max_perturbed_pivots=',max_perturbed_pivots
    passed=phi_l2(CONVERGENCE_LEVELS) < 2.0e-5_dp .and. &
      phi_h1(CONVERGENCE_LEVELS) < 2.0e-3_dp .and. &
      u_l2(CONVERGENCE_LEVELS) < 2.0e-5_dp .and. &
      u_h1(CONVERGENCE_LEVELS) < 2.0e-3_dp .and. &
      phi_l2_order > 2.7_dp .and. phi_h1_order > 1.8_dp .and. &
      u_l2_order > 2.7_dp .and. u_h1_order > 1.8_dp .and. max_perturbed_pivots == 0
  end subroutine check_v02

  subroutine check_v03(passed,status)
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(M0CoefficientsType) :: coefficients
    type(MeshType) :: mesh
    type(PardisoDiagnosticsType) :: diagnostics
    real(dp), allocatable :: solution(:)
    real(dp) :: eta_l2,eta_h1,current_error,right_min,tolerance
    integer :: node

    call clear_status(status)
    call solve_m0_benchmark(4,M0_V03,coefficients,mesh,solution,diagnostics,status)
    if (.not. status%is_ok()) return
    call compute_m0_field_error(mesh,solution,M0_V03,F_ETA_N,coefficients, &
      eta_l2,eta_h1,status)
    if (.not. status%is_ok()) return
    call compute_m0_max_current(mesh,solution,M0_V03,coefficients,current_error,status)
    if (.not. status%is_ok()) return
    tolerance=1.0e-10_dp*coefficients%length
    right_min=huge(1.0_dp)
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(1,node)-coefficients%length) <= tolerance) then
        right_min=min(right_min,solution(global_dof(node,F_ETA_N)))
      end if
    end do
    call print_metric('v03_eta_n_l2_error',eta_l2)
    call print_metric('v03_eta_n_h1_error',eta_h1)
    call print_metric('v03_normalized_max_current',current_error)
    call print_metric('v03_right_face_min_eta_n',right_min)
    write(*,'(a,i0)') 'v03_perturbed_pivots=',diagnostics%perturbed_pivots
    passed=eta_l2 < 1.0e-10_dp .and. eta_h1 < 1.0e-10_dp .and. &
      current_error < 1.0e-10_dp .and. right_min > 0.0_dp .and. &
      diagnostics%perturbed_pivots == 0
  end subroutine check_v03

  subroutine check_v04(passed,status)
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(M0CoefficientsType) :: coefficients
    type(MeshType) :: mesh
    type(PardisoDiagnosticsType) :: diagnostics
    real(dp), allocatable :: solution(:)
    real(dp) :: eta_l2,eta_h1,current_error,right_max,tolerance
    integer :: node

    call clear_status(status)
    call solve_m0_benchmark(4,M0_V04,coefficients,mesh,solution,diagnostics,status)
    if (.not. status%is_ok()) return
    call compute_m0_field_error(mesh,solution,M0_V04,F_ETA_P,coefficients, &
      eta_l2,eta_h1,status)
    if (.not. status%is_ok()) return
    call compute_m0_max_current(mesh,solution,M0_V04,coefficients,current_error,status)
    if (.not. status%is_ok()) return
    tolerance=1.0e-10_dp*coefficients%length
    right_max=-huge(1.0_dp)
    do node=1,mesh%number_of_nodes()
      if (abs(mesh%coordinates(1,node)-coefficients%length) <= tolerance) then
        right_max=max(right_max,solution(global_dof(node,F_ETA_P)))
      end if
    end do
    call print_metric('v04_eta_p_l2_error',eta_l2)
    call print_metric('v04_eta_p_h1_error',eta_h1)
    call print_metric('v04_normalized_max_current',current_error)
    call print_metric('v04_right_face_max_eta_p',right_max)
    write(*,'(a,i0)') 'v04_perturbed_pivots=',diagnostics%perturbed_pivots
    passed=eta_l2 < 1.0e-10_dp .and. eta_h1 < 1.0e-10_dp .and. &
      current_error < 1.0e-10_dp .and. right_max < 0.0_dp .and. &
      diagnostics%perturbed_pivots == 0
  end subroutine check_v04

  subroutine check_v05(passed,status)
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(M0CoefficientsType) :: coefficients
    type(MeshType) :: mesh
    type(PardisoDiagnosticsType) :: diagnostics
    real(dp), allocatable :: solution(:)
    real(dp) :: linear_l2,linear_h1,uniform_l2,uniform_h1,heat_balance
    integer :: maximum_perturbed_pivots

    call clear_status(status)
    maximum_perturbed_pivots=0
    call solve_m0_benchmark(4,M0_V05_LINEAR,coefficients,mesh,solution,diagnostics,status)
    if (.not. status%is_ok()) return
    maximum_perturbed_pivots=max(maximum_perturbed_pivots,diagnostics%perturbed_pivots)
    call compute_m0_field_error(mesh,solution,M0_V05_LINEAR,F_THETA,coefficients, &
      linear_l2,linear_h1,status)
    if (.not. status%is_ok()) return
    call solve_m0_benchmark(4,M0_V05_UNIFORM,coefficients,mesh,solution,diagnostics,status)
    if (.not. status%is_ok()) return
    maximum_perturbed_pivots=max(maximum_perturbed_pivots,diagnostics%perturbed_pivots)
    call compute_m0_field_error(mesh,solution,M0_V05_UNIFORM,F_THETA,coefficients, &
      uniform_l2,uniform_h1,status)
    if (.not. status%is_ok()) return
    call compute_m0_uniform_heat_balance(mesh,solution,coefficients,heat_balance,status)
    if (.not. status%is_ok()) return
    call print_metric('v05_linear_temperature_l2_error',linear_l2)
    call print_metric('v05_linear_temperature_h1_error',linear_h1)
    call print_metric('v05_uniform_source_l2_error',uniform_l2)
    call print_metric('v05_uniform_source_h1_error',uniform_h1)
    call print_metric('v05_uniform_source_energy_balance',heat_balance)
    write(*,'(a,i0)') 'v05_max_perturbed_pivots=',maximum_perturbed_pivots
    passed=max(linear_l2,linear_h1,uniform_l2,uniform_h1) < 1.0e-11_dp .and. &
      heat_balance < 1.0e-11_dp .and. maximum_perturbed_pivots == 0
  end subroutine check_v05

  subroutine check_v06(passed,status)
    logical, intent(out) :: passed
    type(SolverStatusType), intent(out) :: status
    type(M0CoefficientsType) :: coefficients
    type(MeshType) :: mesh
    type(PardisoDiagnosticsType) :: diagnostics
    real(dp), allocatable :: solution(:)
    real(dp) :: phi_l2(CONVERGENCE_LEVELS),phi_h1(CONVERGENCE_LEVELS)
    real(dp) :: u_l2(CONVERGENCE_LEVELS),u_h1(CONVERGENCE_LEVELS)
    real(dp) :: eta_l2(CONVERGENCE_LEVELS),eta_h1(CONVERGENCE_LEVELS)
    real(dp) :: theta_l2(CONVERGENCE_LEVELS),theta_h1(CONVERGENCE_LEVELS)
    real(dp) :: phi_l2_order,phi_h1_order,current_error
    integer :: level,max_perturbed_pivots

    call clear_status(status)
    max_perturbed_pivots=0
    do level=1,CONVERGENCE_LEVELS
      call solve_m0_benchmark(ELEMENT_COUNTS(level),M0_V06,coefficients,mesh,solution, &
        diagnostics,status)
      if (.not. status%is_ok()) return
      max_perturbed_pivots=max(max_perturbed_pivots,diagnostics%perturbed_pivots)
      call compute_m0_field_error(mesh,solution,M0_V06,F_PHI,coefficients, &
        phi_l2(level),phi_h1(level),status)
      if (.not. status%is_ok()) return
      call compute_m0_field_error(mesh,solution,M0_V06,F_U1,coefficients, &
        u_l2(level),u_h1(level),status)
      if (.not. status%is_ok()) return
      call compute_m0_field_error(mesh,solution,M0_V06,F_ETA_N,coefficients, &
        eta_l2(level),eta_h1(level),status)
      if (.not. status%is_ok()) return
      call compute_m0_field_error(mesh,solution,M0_V06,F_THETA,coefficients, &
        theta_l2(level),theta_h1(level),status)
      if (.not. status%is_ok()) return
      call print_level_metrics('v06',ELEMENT_COUNTS(level),phi_l2(level),phi_h1(level))
    end do
    call compute_m0_max_current(mesh,solution,M0_V06,coefficients,current_error,status)
    if (.not. status%is_ok()) return
    phi_l2_order=observed_order(phi_l2)
    phi_h1_order=observed_order(phi_h1)
    call print_metric('v06_phi_l2_order',phi_l2_order)
    call print_metric('v06_phi_h1_order',phi_h1_order)
    call print_metric('v06_u_l2_final',u_l2(CONVERGENCE_LEVELS))
    call print_metric('v06_u_h1_final',u_h1(CONVERGENCE_LEVELS))
    call print_metric('v06_eta_n_l2_final',eta_l2(CONVERGENCE_LEVELS))
    call print_metric('v06_eta_n_h1_final',eta_h1(CONVERGENCE_LEVELS))
    call print_metric('v06_theta_l2_final',theta_l2(CONVERGENCE_LEVELS))
    call print_metric('v06_theta_h1_final',theta_h1(CONVERGENCE_LEVELS))
    call print_metric('v06_normalized_max_current',current_error)
    write(*,'(a,i0)') 'v06_max_perturbed_pivots=',max_perturbed_pivots
    passed=phi_l2(CONVERGENCE_LEVELS) < 2.0e-5_dp .and. &
      phi_h1(CONVERGENCE_LEVELS) < 2.0e-3_dp .and. &
      max(u_l2(CONVERGENCE_LEVELS),eta_l2(CONVERGENCE_LEVELS), &
      theta_l2(CONVERGENCE_LEVELS)) < 2.0e-5_dp .and. &
      max(u_h1(CONVERGENCE_LEVELS),eta_h1(CONVERGENCE_LEVELS), &
      theta_h1(CONVERGENCE_LEVELS)) < 2.0e-3_dp .and. &
      phi_l2_order > 2.7_dp .and. phi_h1_order > 1.8_dp .and. &
      current_error < 1.0e-8_dp .and. max_perturbed_pivots == 0
  end subroutine check_v06

  real(dp) function observed_order(errors) result(order)
    real(dp), intent(in) :: errors(:)
    integer :: count
    count=size(errors)
    if (count < 2 .or. errors(count) <= 0.0_dp .or. errors(count-1) <= 0.0_dp) then
      order=-huge(1.0_dp)
    else
      order=log(errors(count-1)/errors(count))/log(2.0_dp)
    end if
  end function observed_order

  subroutine print_level_metrics(prefix,elements,l2_error,h1_error)
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: elements
    real(dp), intent(in) :: l2_error,h1_error
    character(len=64) :: label
    write(label,'(a,"_ne",i0,"_phi_l2")') trim(prefix),elements
    call print_metric(trim(label),l2_error)
    write(label,'(a,"_ne",i0,"_phi_h1")') trim(prefix),elements
    call print_metric(trim(label),h1_error)
  end subroutine print_level_metrics

  subroutine print_metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine print_metric

  subroutine copy_failure(source,destination)
    type(SolverStatusType), intent(in) :: source
    type(SolverStatusType), intent(out) :: destination
    call set_status(destination,source%code,trim(source%message))
  end subroutine copy_failure

  pure function uppercase(text) result(converted)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: converted
    integer :: index,code
    do index=1,len(text)
      code=iachar(text(index:index))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        converted(index:index)=achar(code-32)
      else
        converted(index:index)=text(index:index)
      end if
    end do
  end function uppercase

end module analytic_benchmarks_m

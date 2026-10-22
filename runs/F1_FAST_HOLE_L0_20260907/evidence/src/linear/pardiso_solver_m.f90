module pardiso_solver_m
  use kinds_m, only : dp, i8
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INVALID_ARGUMENT, STATUS_LINEAR_SOLVER_FAILED
  implicit none
  private

  integer, parameter :: PARDISO_REAL_NONSYMMETRIC = 11

  type, public :: PardisoDiagnosticsType
    integer :: error_code = 0
    integer :: factor_nonzeros = 0
    integer :: perturbed_pivots = 0
    integer :: iterative_refinements = 0
  end type PardisoDiagnosticsType

  public :: solve_csr_pardiso

  interface
    subroutine pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs, &
                       iparm,msglvl,b,x,error)
      import :: dp, i8
      integer(i8), intent(inout) :: pt(*)
      integer, intent(in) :: maxfct,mnum,mtype,phase,n,nrhs,msglvl
      integer, intent(in) :: ia(*),ja(*)
      integer, intent(inout) :: perm(*),iparm(*)
      real(dp), intent(in) :: a(*)
      real(dp), intent(inout) :: b(*)
      real(dp), intent(out) :: x(*)
      integer, intent(out) :: error
    end subroutine pardiso
  end interface

contains

  subroutine solve_csr_pardiso(number_of_equations,ia,ja,values,rhs,solution, &
                               diagnostics,status)
    integer, intent(in) :: number_of_equations
    integer, intent(in) :: ia(:),ja(:)
    real(dp), intent(in) :: values(:),rhs(:)
    real(dp), intent(out) :: solution(:)
    type(PardisoDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    integer(i8) :: pt(64)
    integer :: iparm(64),maxfct,mnum,mtype,phase,nrhs,msglvl,error,release_error
    integer, allocatable :: perm(:)
    real(dp), allocatable :: rhs_work(:)
    character(len=16) :: storage_mode
    integer :: environment_status

    call clear_status(status)
    diagnostics = PardisoDiagnosticsType()
    solution = 0.0_dp
    if (number_of_equations < 1 .or. size(ia) /= number_of_equations+1 .or. &
        size(rhs) /= number_of_equations .or. size(solution) /= number_of_equations) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'PARDISO array size mismatch.')
      return
    end if
    if (ia(1) /= 1 .or. ia(number_of_equations+1)-1 /= size(values) .or. &
        size(ja) /= size(values)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'PARDISO requires consistent one-based CSR arrays.')
      return
    end if
    if (any(ia(2:) < ia(:number_of_equations)) .or. any(ja < 1) .or. &
        any(ja > number_of_equations)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'PARDISO CSR indices are invalid.')
      return
    end if

    allocate(perm(number_of_equations),rhs_work(number_of_equations))
    pt = 0_i8
    iparm = 0
    perm = 0
    rhs_work = rhs
    maxfct = 1
    mnum = 1
    mtype = PARDISO_REAL_NONSYMMETRIC
    nrhs = 1
    msglvl = 0

    ! Freeze explicit, deterministic LP64 settings instead of relying on MKL defaults.
    iparm(1) = 1
    iparm(2) = 2
    iparm(8) = 2
    iparm(10) = 13
    iparm(11) = 1
    iparm(13) = 1
    iparm(18) = -1
    iparm(27) = 1
    iparm(35) = 0

    ! Optional bounded-memory storage only; same sparse equations and tolerances.
    ! Fortran iparm(60)=2 is the documented oneMKL out-of-core mode.
    storage_mode='IC'
    call get_environment_variable('TPS_PARDISO_STORAGE',storage_mode,status=environment_status)
    if (environment_status==1) storage_mode='IC'
    if (trim(storage_mode)=='OOC') then
      iparm(60)=2
    else if (trim(storage_mode)/='IC') then
      call set_status(status,STATUS_INVALID_ARGUMENT,'TPS_PARDISO_STORAGE must be IC or OOC.');return
    end if

    phase = 13
    call pardiso(pt,maxfct,mnum,mtype,phase,number_of_equations,values,ia,ja, &
      perm,nrhs,iparm,msglvl,rhs_work,solution,error)
    diagnostics%error_code = error
    diagnostics%iterative_refinements = iparm(7)
    diagnostics%perturbed_pivots = iparm(14)
    diagnostics%factor_nonzeros = iparm(18)

    phase = -1
    release_error = 0
    call pardiso(pt,maxfct,mnum,mtype,phase,number_of_equations,values,ia,ja, &
      perm,nrhs,iparm,msglvl,rhs_work,solution,release_error)

    if (error /= 0) then
      call set_status(status,STATUS_LINEAR_SOLVER_FAILED, &
        'PARDISO phase 13 failed with error '//integer_to_string(error)//'.')
      return
    end if
    if (release_error /= 0) then
      call set_status(status,STATUS_LINEAR_SOLVER_FAILED, &
        'PARDISO release failed with error '//integer_to_string(release_error)//'.')
    end if
  end subroutine solve_csr_pardiso

  pure function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

end module pardiso_solver_m

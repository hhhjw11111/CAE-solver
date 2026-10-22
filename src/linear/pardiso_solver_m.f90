module pardiso_solver_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : output_unit
  use kinds_m, only : dp, i8
  use field_layout_m, only : NDOF_NODE,F_ETA_P
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
    logical :: split_attempted = .false.
    logical :: split_accepted = .false.
    real(dp) :: split_backward_error = huge(1.0_dp)
    integer :: split_primary_entries = 0
    integer :: split_hole_entries = 0
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
                               diagnostics,status,try_hole_split)
    integer, intent(in) :: number_of_equations
    integer, intent(in) :: ia(:),ja(:)
    real(dp), intent(in) :: values(:),rhs(:)
    real(dp), intent(out) :: solution(:)
    type(PardisoDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    logical, intent(in), optional :: try_hole_split
    type(PardisoDiagnosticsType) :: split_diagnostics
    character(len=16) :: split_mode
    integer :: environment_status
    logical :: enabled,eligible

    split_mode='0'
    call get_environment_variable('TPS_PARDISO_SPLIT_HOLES',split_mode,status=environment_status)
    enabled=environment_status==0.and.trim(split_mode)=='1'
    if(present(try_hole_split)) enabled=try_hole_split
    eligible=.false.
    if(enabled.and.number_of_equations>=2*NDOF_NODE) then
      if(mod(number_of_equations,NDOF_NODE)==0.and.size(ia)==number_of_equations+1.and. &
         size(rhs)==number_of_equations.and.size(solution)==number_of_equations.and.size(ja)==size(values)) then
        if(ia(1)==1.and.ia(number_of_equations+1)==size(values)+1) then
          eligible=all(ia(2:)>=ia(:number_of_equations)).and.all(ja>=1).and.all(ja<=number_of_equations)
        end if
      end if
    end if
    if(eligible) then
      call attempt_hole_split(number_of_equations,ia,ja,values,rhs,solution,split_diagnostics,status)
      if(split_diagnostics%split_accepted) then
        diagnostics=split_diagnostics
        if(environment_status==0.and.trim(split_mode)=='1') then
          write(output_unit,'(a,es12.4,a)') 'pardiso_hole_split=ACCEPTED,full_backward_error=', &
            diagnostics%split_backward_error,','
          write(output_unit,'(a,i0,a,i0)') 'pardiso_primary_entries=',diagnostics%split_primary_entries, &
            ',hole_entries=',diagnostics%split_hole_entries
          flush(output_unit)
        end if
        return
      end if
    end if
    call solve_csr_pardiso_plain(number_of_equations,ia,ja,values,rhs,solution,diagnostics,status)
    if(eligible) then
      diagnostics%split_attempted=.true.
      diagnostics%split_backward_error=split_diagnostics%split_backward_error
      if(environment_status==0.and.trim(split_mode)=='1') then
        write(output_unit,'(a,es12.4,a)') 'pardiso_hole_split=FULL_MATRIX_FALLBACK,rejected_backward_error=', &
          diagnostics%split_backward_error,','
        flush(output_unit)
      end if
    end if
  end subroutine solve_csr_pardiso

  subroutine attempt_hole_split(n,ia,ja,values,rhs,solution,diagnostics,status)
    ! One block-triangular candidate [A 0; C D]. The neglected B block is
    ! checked using every row of the ORIGINAL matrix, including constraints.
    ! This is a linear-solver optimization, never a change to the nonlinear PDE.
    integer, intent(in) :: n,ia(:),ja(:)
    real(dp), intent(in) :: values(:),rhs(:)
    real(dp), intent(out) :: solution(:)
    type(PardisoDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    integer, allocatable :: primary(:),holes(:),ia_p(:),ja_p(:),ia_h(:),ja_h(:)
    real(dp), allocatable :: a_p(:),a_h(:),b_p(:),b_h(:),x_p(:),x_h(:)
    type(PardisoDiagnosticsType) :: diag_p,diag_h
    integer :: row,col,k,np,nh,kp,kh
    real(dp) :: feedback,row_sum,action,denominator,backward
    call clear_status(status)
    diagnostics=PardisoDiagnosticsType();diagnostics%split_attempted=.true.
    solution=0
    if(.not.all(ieee_is_finite(values)).or..not.all(ieee_is_finite(rhs))) return
    allocate(primary(n),holes(n));np=0;nh=0;primary=0;holes=0
    do row=1,n
      if(mod(row-1,NDOF_NODE)+1==F_ETA_P) then
        nh=nh+1;holes(row)=nh
      else
        np=np+1;primary(row)=np
      end if
    end do
    kp=0;kh=0
    do row=1,n
      feedback=0;row_sum=0
      do k=ia(row),ia(row+1)-1
        col=ja(k)
        if(primary(row)>0) then
          row_sum=row_sum+abs(values(k))
          if(primary(col)>0) then
            if(values(k)/=0.0_dp.or.col==row) kp=kp+1
          else
            feedback=feedback+abs(values(k))
          end if
        else if(holes(col)>0) then
          if(values(k)/=0.0_dp.or.col==row) kh=kh+1
        end if
      end do
      ! Avoid two unnecessary factorizations when hole feedback is substantial.
      if(primary(row)>0.and.feedback>1.0e-10_dp*max(row_sum,tiny(1.0_dp))) return
    end do
    ! Remove only exact off-diagonal zeros from the two candidate CSR graphs.
    ! Keep every diagonal, even zero. No tolerance-based coefficient dropping;
    ! the complete original matrix still checks the candidate after both solves.
    diagnostics%split_primary_entries=kp;diagnostics%split_hole_entries=kh
    allocate(ia_p(np+1),ja_p(kp),a_p(kp),b_p(np),x_p(np), &
             ia_h(nh+1),ja_h(kh),a_h(kh),b_h(nh),x_h(nh))
    kp=0;kh=0
    do row=1,n
      if(primary(row)>0) then
        ia_p(primary(row))=kp+1;b_p(primary(row))=rhs(row)
        do k=ia(row),ia(row+1)-1
          col=ja(k)
          if(primary(col)==0) cycle
          if(values(k)==0.0_dp.and.col/=row) cycle
          kp=kp+1;ja_p(kp)=primary(col);a_p(kp)=values(k)
        end do
      else
        ia_h(holes(row))=kh+1
        do k=ia(row),ia(row+1)-1
          col=ja(k)
          if(holes(col)==0) cycle
          if(values(k)==0.0_dp.and.col/=row) cycle
          kh=kh+1;ja_h(kh)=holes(col);a_h(kh)=values(k)
        end do
      end if
    end do
    ia_p(np+1)=kp+1;ia_h(nh+1)=kh+1
    call solve_csr_pardiso_plain(np,ia_p,ja_p,a_p,b_p,x_p,diag_p,status)
    if(.not.status%is_ok().or..not.all(ieee_is_finite(x_p))) return
    do row=1,n
      if(primary(row)>0) then
        solution(row)=x_p(primary(row))
      else
        b_h(holes(row))=rhs(row)
        do k=ia(row),ia(row+1)-1
          col=ja(k)
          if(primary(col)>0) b_h(holes(row))=b_h(holes(row))-values(k)*x_p(primary(col))
        end do
      end if
    end do
    call solve_csr_pardiso_plain(nh,ia_h,ja_h,a_h,b_h,x_h,diag_h,status)
    if(.not.status%is_ok().or..not.all(ieee_is_finite(x_h))) return
    do row=1,n
      if(holes(row)>0) solution(row)=x_h(holes(row))
    end do
    ! Recover exact single-diagonal constraints before the componentwise check:
    ! for a zero RHS, any roundoff in x would otherwise give relative error 1.
    do row=1,n
      kp=0;kh=0
      do k=ia(row),ia(row+1)-1
        if(values(k)==0) cycle
        kp=kp+1;kh=k
      end do
      if(kp==1) then
        if(ja(kh)==row) solution(row)=rhs(row)/values(kh)
      end if
    end do
    backward=0
    do row=1,n
      action=0;denominator=abs(rhs(row))
      do k=ia(row),ia(row+1)-1
        action=action+values(k)*solution(ja(k))
        denominator=denominator+abs(values(k)*solution(ja(k)))
      end do
      if(.not.ieee_is_finite(action).or..not.ieee_is_finite(denominator)) return
      backward=max(backward,abs(action-rhs(row))/max(denominator,tiny(1.0_dp)))
    end do
    diagnostics%split_backward_error=backward
    if(.not.ieee_is_finite(backward).or.backward>1.0e-11_dp) return
    diagnostics%split_accepted=.true.
    diagnostics%factor_nonzeros=diag_p%factor_nonzeros+diag_h%factor_nonzeros
    diagnostics%perturbed_pivots=diag_p%perturbed_pivots+diag_h%perturbed_pivots
    diagnostics%iterative_refinements=max(diag_p%iterative_refinements,diag_h%iterative_refinements)
  end subroutine attempt_hole_split

  subroutine solve_csr_pardiso_plain(number_of_equations,ia,ja,values,rhs,solution, &
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
  end subroutine solve_csr_pardiso_plain

  pure function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

end module pardiso_solver_m

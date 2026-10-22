program pardiso_representation_driver
  ! Diagnostic only. No PDE or production solver source is changed.
  ! The direct call settings below are copied from the frozen production wrapper.
  ! The linked ORIGINAL wrapper is called separately to verify the copy.
  use kinds_m, only : dp,i8
  use status_m, only : SolverStatusType
  use pardiso_solver_m, only : solve_csr_pardiso,PardisoDiagnosticsType
  implicit none
  integer(i8) :: pt(64)
  integer :: iparm(64),before(64),after(64),n,nz,u,error,release_error
  integer :: maxfct,mnum,mtype,phase,nrhs,msglvl
  integer,allocatable :: ia(:),ja(:),perm(:)
  real(dp),allocatable :: a(:),rhs(:),work(:),x(:),saved_x(:),wrapper_x(:)
  character(len=2048) :: input_path,output_path
  character(len=16) :: storage_mode
  type(SolverStatusType) :: status
  type(PardisoDiagnosticsType) :: diagnostics
  interface
    subroutine pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,b,x,error)
      import :: dp,i8
      integer(i8),intent(inout) :: pt(*)
      integer,intent(in) :: maxfct,mnum,mtype,phase,n,nrhs,msglvl
      integer,intent(in) :: ia(*),ja(*)
      integer,intent(inout) :: perm(*),iparm(*)
      real(dp),intent(in) :: a(*)
      real(dp),intent(inout) :: b(*)
      real(dp),intent(out) :: x(*)
      integer,intent(out) :: error
    end subroutine pardiso
  end interface
  call get_command_argument(1,input_path)
  call get_command_argument(2,output_path)
  storage_mode='IC'
  call get_environment_variable('TPS_PARDISO_STORAGE',storage_mode)
  if(trim(storage_mode)/='IC') error stop 'This diagnostic freezes IC, matching matrix export.'
  open(newunit=u,file=trim(input_path),access='stream',form='unformatted',status='old',action='read')
  read(u) n,nz
  allocate(ia(n+1),ja(nz),a(nz),rhs(n),work(n),x(n),saved_x(n),wrapper_x(n),perm(n))
  read(u) ia,ja,a,rhs
  close(u)
  pt=0_i8;iparm=0;perm=0;work=rhs;x=0.0_dp
  maxfct=1;mnum=1;mtype=11;nrhs=1;msglvl=0
  iparm(1)=1
  iparm(2)=2
  iparm(8)=2
  iparm(10)=13
  iparm(11)=1
  iparm(13)=1
  iparm(18)=-1
  iparm(27)=1
  iparm(35)=0
  before=iparm
  phase=13
  call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,work,x,error)
  after=iparm;saved_x=x
  phase=-1;release_error=0
  call pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,work,x,release_error)
  call solve_csr_pardiso(n,ia,ja,a,rhs,wrapper_x,diagnostics,status,try_hole_split=.false.)
  open(newunit=u,file=trim(output_path),access='stream',form='unformatted',status='new',action='write')
  write(u) 1,n,nz,error,release_error,status%code,diagnostics%error_code, &
    diagnostics%perturbed_pivots,diagnostics%iterative_refinements,diagnostics%factor_nonzeros
  write(u) before,after,iparm,saved_x,wrapper_x
  close(u)
  write(*,'(a,i0,a,i0,a,i0,a,i0)') 'n=',n,',nnz=',nz,',error=',error,',pivots=',after(14)
  write(*,'(a,i0,a,es24.16)') 'wrapper_pivots=',diagnostics%perturbed_pivots, &
    ',direct_wrapper_maxdiff=',maxval(abs(saved_x-wrapper_x))
  if(error/=0.or.release_error/=0.or..not.status%is_ok()) stop 1
end program pardiso_representation_driver

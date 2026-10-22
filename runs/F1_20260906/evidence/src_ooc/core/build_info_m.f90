module build_info_m
  implicit none
  private

  character(len=*), parameter, public :: SOLVER_NAME = 'TPSolverV2'
  character(len=*), parameter, public :: SOLVER_VERSION = '0.14.1-f1-bounded-memory'
  character(len=*), parameter, public :: THEORY_VERSION = '2.1'
  character(len=*), parameter, public :: ARCHITECTURE_ID = 'TPS-V2-ARCH-001'
  character(len=*), parameter, public :: VERIFIED_SCOPE = 'See run-specific evidence; legacy graded-hole geometry INVALID'

  public :: print_build_info

contains

  subroutine print_build_info(unit_id)
    integer, intent(in) :: unit_id
    write(unit_id,'(a)') 'solver=' // SOLVER_NAME
    write(unit_id,'(a)') 'version=' // SOLVER_VERSION
    write(unit_id,'(a)') 'theory_version=' // THEORY_VERSION
    write(unit_id,'(a)') 'architecture_id=' // ARCHITECTURE_ID
    write(unit_id,'(a)') 'verified_scope=' // VERIFIED_SCOPE
    write(unit_id,'(a)') 'dof_order=u1,u2,u3,phi,eta_n,eta_p,theta'
  end subroutine print_build_info

end module build_info_m

module status_m
  implicit none
  private

  integer, parameter, public :: STATUS_OK = 0
  integer, parameter, public :: STATUS_INVALID_ARGUMENT = 10
  integer, parameter, public :: STATUS_INPUT_ERROR = 20
  integer, parameter, public :: STATUS_GEOMETRY_ERROR = 30
  integer, parameter, public :: STATUS_NOT_IMPLEMENTED = 40
  integer, parameter, public :: STATUS_VERIFICATION_FAILED = 50
  integer, parameter, public :: STATUS_LINEAR_SOLVER_FAILED = 60
  integer, parameter, public :: STATUS_NONLINEAR_SOLVER_FAILED = 70
  integer, parameter, public :: STATUS_CONTINUATION_FAILED = 80

  type, public :: SolverStatusType
    integer :: code = STATUS_OK
    character(len=512) :: message = 'OK'
  contains
    procedure :: is_ok => status_is_ok
  end type SolverStatusType

  public :: clear_status, set_status

contains

  logical function status_is_ok(this)
    class(SolverStatusType), intent(in) :: this
    status_is_ok = this%code == STATUS_OK
  end function status_is_ok

  subroutine clear_status(status)
    type(SolverStatusType), intent(out) :: status
    status%code = STATUS_OK
    status%message = 'OK'
  end subroutine clear_status

  subroutine set_status(status, code, message)
    type(SolverStatusType), intent(out) :: status
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    status%code = code
    status%message = ''
    status%message = trim(message)
  end subroutine set_status

end module status_m

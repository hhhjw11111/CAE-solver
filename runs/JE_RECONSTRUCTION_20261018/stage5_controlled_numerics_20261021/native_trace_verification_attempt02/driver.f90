program trace_controlled_verification
  use status_m, only: SolverStatusType
  use supg_trace_iso_check_m, only: verify_supg_trace_iso
  use supg_trace_thin_check_m, only: verify_supg_trace_thin
  implicit none
  type(SolverStatusType) :: status
  call verify_supg_trace_iso(status)
  if (.not.status%is_ok()) then
    write(*,'(a)') trim(status%message)
    stop 1
  end if
  call verify_supg_trace_thin(status)
  if (.not.status%is_ok()) then
    write(*,'(a)') trim(status%message)
    stop 2
  end if
  write(*,'(a)') 'trace_controlled_verification=PASS'
end program trace_controlled_verification

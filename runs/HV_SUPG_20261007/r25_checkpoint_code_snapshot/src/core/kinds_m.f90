module kinds_m
  use, intrinsic :: iso_fortran_env, only : real64, int32, int64
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: i4 = int32
  integer, parameter, public :: i8 = int64
end module kinds_m

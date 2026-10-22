module constants_m
  use kinds_m, only : dp
  implicit none
  private

  real(dp), parameter, public :: ELEMENTARY_CHARGE = 1.602176634e-19_dp
  real(dp), parameter, public :: BOLTZMANN_CONSTANT = 1.380649e-23_dp
  real(dp), parameter, public :: ABSOLUTE_ZERO_TOL = 1.0e-12_dp
  real(dp), parameter, public :: GEOMETRY_TOL = 1.0e-12_dp
  real(dp), parameter, public :: VERIFICATION_TOL = 5.0e-12_dp
end module constants_m

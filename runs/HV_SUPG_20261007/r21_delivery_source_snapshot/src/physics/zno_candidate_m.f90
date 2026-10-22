module zno_candidate_m
  use kinds_m, only : dp
  use model_types_m, only : MaterialType
  implicit none
  private
  public :: define_zno_candidate
contains
  subroutine define_zno_candidate(model,material)
    character(len=*), intent(in) :: model
    type(MaterialType), intent(out) :: material
    real(dp) :: alpha(6)
    integer :: i
    material=MaterialType()
    material%name='ZNO_300K_CANDIDATE_V2_'//trim(model)
    ! Candidate effective inputs; mixed-source, not same-specimen calibration.
    ! c||+z, Voigt xx,yy,zz,yz,zx,xy with engineering shear strain.
    material%elastic(1,1)=209.7e9_dp;material%elastic(2,2)=209.7e9_dp
    material%elastic(3,3)=210.9e9_dp
    material%elastic(1,2)=121.1e9_dp;material%elastic(2,1)=121.1e9_dp
    material%elastic(1,3)=105.1e9_dp;material%elastic(3,1)=105.1e9_dp
    material%elastic(2,3)=105.1e9_dp;material%elastic(3,2)=105.1e9_dp
    material%elastic(4,4)=42.47e9_dp;material%elastic(5,5)=42.47e9_dp
    material%elastic(6,6)=44.3e9_dp
    material%piezoelectric(1,5)=-0.37_dp;material%piezoelectric(2,4)=-0.37_dp
    material%piezoelectric(3,1:3)=[-0.62_dp,-0.62_dp,0.96_dp]
    do i=1,3
      material%permittivity(i,i)=8.54_dp*8.8541878128e-12_dp
      material%thermal_conductivity(i,i)=60.0_dp
      material%mobility_n0(i,i)=0.0205_dp
      material%mobility_p0(i,i)=0.005_dp
    end do
    material%permittivity(3,3)=10.2_dp*8.8541878128e-12_dp
    alpha=[4.31e-6_dp,4.31e-6_dp,2.49e-6_dp,0.0_dp,0.0_dp,0.0_dp]
    material%thermal_stress=matmul(material%elastic,alpha)
    material%pyroelectric=[0.0_dp,0.0_dp,-6.9e-6_dp]
    material%beta_n=2.0_dp;material%beta_p=2.0_dp
    material%reference_temperature=300.0_dp
    material%reference_n=6.0e22_dp
    material%intrinsic_concentration=2.9e-4_dp
    material%reference_p=material%intrinsic_concentration**2/material%reference_n
    material%donor_concentration=material%reference_n-material%reference_p
    material%frozen_transport=trim(model)=='N0' .or. trim(model)=='LQ'
    material%linear_carriers=trim(model)=='LQ'
  end subroutine define_zno_candidate
end module zno_candidate_m

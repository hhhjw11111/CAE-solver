module reference_material_m
  use kinds_m, only : dp
  use model_types_m, only : MaterialType
  implicit none
  private

  public :: define_tps_g08_reference_material

contains

  subroutine define_tps_g08_reference_material(material)
    type(MaterialType), intent(out) :: material

    material=MaterialType()
    material%name='TPS_G08_REFERENCE'
    material%elastic=0.0_dp
    material%elastic(1,1)=1.20e11_dp
    material%elastic(2,2)=1.10e11_dp
    material%elastic(3,3)=1.35e11_dp
    material%elastic(4,4)=4.20e10_dp
    material%elastic(5,5)=4.60e10_dp
    material%elastic(6,6)=4.00e10_dp
    material%elastic(1,2)=4.0e10_dp; material%elastic(2,1)=4.0e10_dp
    material%elastic(1,3)=3.5e10_dp; material%elastic(3,1)=3.5e10_dp
    material%elastic(2,3)=3.8e10_dp; material%elastic(3,2)=3.8e10_dp
    material%piezoelectric=reshape([ &
      0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.15_dp,-0.12_dp, &
      0.0_dp,0.0_dp,0.0_dp,0.14_dp,0.0_dp,0.10_dp, &
      -0.30_dp,-0.28_dp,0.65_dp,0.0_dp,0.0_dp,0.0_dp ],[3,6])
    material%permittivity=0.0_dp
    material%permittivity(1,1)=8.0e-11_dp
    material%permittivity(2,2)=8.5e-11_dp
    material%permittivity(3,3)=9.2e-11_dp
    material%thermal_stress=[1.8e6_dp,1.7e6_dp,2.1e6_dp,0.0_dp,0.0_dp,0.0_dp]
    material%pyroelectric=[0.2e-5_dp,-0.1e-5_dp,1.5e-5_dp]
    material%thermal_conductivity=0.0_dp
    material%thermal_conductivity(1,1)=18.0_dp
    material%thermal_conductivity(2,2)=20.0_dp
    material%thermal_conductivity(3,3)=24.0_dp
    material%mobility_n0=0.0_dp
    material%mobility_n0(1,1)=0.020_dp
    material%mobility_n0(2,2)=0.024_dp
    material%mobility_n0(3,3)=0.028_dp
    material%mobility_p0=0.0_dp
    material%mobility_p0(1,1)=0.006_dp
    material%mobility_p0(2,2)=0.007_dp
    material%mobility_p0(3,3)=0.008_dp
    material%beta_n=1.35_dp
    material%beta_p=1.10_dp
    material%reference_temperature=300.0_dp
    material%intrinsic_concentration=1.0e16_dp
    material%donor_concentration=1.0e21_dp
    material%acceptor_concentration=2.0e20_dp
    material%reference_n=8.0e20_dp
    material%reference_p=1.25e11_dp
  end subroutine define_tps_g08_reference_material

end module reference_material_m

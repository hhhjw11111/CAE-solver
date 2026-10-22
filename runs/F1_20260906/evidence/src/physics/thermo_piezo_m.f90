module thermo_piezo_m
  use kinds_m, only : dp
  use model_types_m, only : MaterialType
  implicit none
  private

  public :: evaluate_thermo_piezo

contains

  pure subroutine evaluate_thermo_piezo(material, strain, electric_field, temperature, &
                                        stress, electric_displacement)
    type(MaterialType), intent(in) :: material
    real(dp), intent(in) :: strain(6), electric_field(3), temperature
    real(dp), intent(out) :: stress(6), electric_displacement(3)
    real(dp) :: temperature_change

    temperature_change = temperature - material%reference_temperature
    stress = matmul(material%elastic, strain) &
           - matmul(transpose(material%piezoelectric), electric_field) &
           - material%thermal_stress * temperature_change
    electric_displacement = matmul(material%piezoelectric, strain) &
                          + matmul(material%permittivity, electric_field) &
                          + material%pyroelectric * temperature_change
  end subroutine evaluate_thermo_piezo

end module thermo_piezo_m

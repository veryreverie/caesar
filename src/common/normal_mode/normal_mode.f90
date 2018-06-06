! ======================================================================
! Types and methods for working with normal mode co-ordinates. Includes:
!    - ComplexMode and RealMode, which store their cartesian displacements.
!    - Displacements, which store displacements in normal mode co-ordinates.
!    - Univariates, Monomials and Polynomials, which store basis functions.
!    - Conversions between complex and real co-ordinates.
! ======================================================================
module normal_mode_module
  use cartesian_displacement_submodule
  use cartesian_force_submodule
  use complex_mode_submodule
  use real_mode_submodule
  use complex_single_mode_displacement_submodule
  use real_single_mode_displacement_submodule
  use complex_mode_displacement_submodule
  use real_mode_displacement_submodule
  use complex_polynomial_submodule
  use real_polynomial_submodule
  use real_complex_conversion_submodule
  use complex_mode_symmetry_submodule
  implicit none
end module

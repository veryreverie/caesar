! ======================================================================
! Types and methods for working with normal mode co-ordinates. Includes:
!    - ComplexMode and RealMode, which store their cartesian displacements.
!    - Displacements, which store displacements in normal mode co-ordinates.
!    - Univariates, Monomials and Polynomials, which store basis functions.
!    - Conversions between complex and real co-ordinates.
! ======================================================================
! This module is simply an interface for the various subsidiary
!    normal mode modules.
module caesar_normal_mode_module
  use caesar_cartesian_displacement_module
  use caesar_cartesian_force_module
  use caesar_cartesian_hessian_module
  use caesar_mass_weighted_displacement_module
  use caesar_mass_weighted_force_module
  use caesar_real_mode_module
  use caesar_complex_mode_module
  use caesar_real_single_mode_displacement_module
  use caesar_real_single_mode_force_module
  use caesar_complex_single_mode_displacement_module
  use caesar_complex_single_mode_force_module
  use caesar_real_mode_displacement_module
  use caesar_real_mode_force_module
  use caesar_complex_mode_displacement_module
  use caesar_complex_mode_force_module
  use caesar_complex_polynomial_module
  use caesar_real_polynomial_module
  use caesar_paired_polynomial_module
  use caesar_real_complex_conversion_module
  use caesar_complex_mode_symmetry_module
  use caesar_degenerate_subspace_module
  use caesar_process_modes_module
  implicit none
end module

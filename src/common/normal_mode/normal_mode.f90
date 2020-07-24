! ======================================================================
! Types and methods for working with normal mode co-ordinates. Includes:
!    - ComplexMode and RealMode, which store their cartesian displacements.
!    - Displacements, which store displacements in normal mode co-ordinates.
!    - Univariates, Monomials and Polynomials, which store basis functions.
!    - Conversions between complex and real co-ordinates.
! ======================================================================
! This module is simply an interface for the various subsidiary
!    normal mode modules.
module normal_mode_module
  use cartesian_displacement_module
  use cartesian_force_module
  use cartesian_hessian_module
  use mass_weighted_displacement_module
  use mass_weighted_force_module
  use real_mode_module
  use complex_mode_module
  use real_single_mode_displacement_module
  use real_single_mode_force_module
  use complex_single_mode_displacement_module
  use complex_single_mode_force_module
  use real_mode_displacement_module
  use real_mode_force_module
  use complex_mode_displacement_module
  use complex_mode_force_module
  use complex_polynomial_module
  use real_polynomial_module
  use paired_polynomial_module
  use real_complex_conversion_module
  use complex_mode_symmetry_module
  use degenerate_subspace_module
  use process_modes_module
  implicit none
end module

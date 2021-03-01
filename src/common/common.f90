! ======================================================================
! Provides higher-level utilities, specific to the function of the code.
! Utilities include:
!    - The Structure type and its subsidiaries, which handle atomic structure
!         information, such as lattice parameters and atomic informaion.
!    - Interfaces to electronic structure codes, including CASTEP and Quip.
!    - Types for handling normal mode co-ordinates.
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module caesar_common_module
  use caesar_utils_module
  
  use caesar_structure_module
  use caesar_normal_mode_module
  use caesar_dynamical_matrices_module
  use caesar_harmonic_stress_module
  use caesar_electronic_structure_module
  use caesar_observables_module
  implicit none
end module

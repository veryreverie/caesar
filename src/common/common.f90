! ======================================================================
! Provides higher-level utilities, specific to the function of the code.
! Utilities include:
!    - The Structure type and its subsidiaries, which handle atomic structure
!         information, such as lattice parameters and atomic informaion.
!    - Interfaces to electronic structure codes, including CASTEP and Quip.
!    - Types for handling normal mode co-ordinates.
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module common_module
  use utils_module
  
  use structure_module
  use electronic_structure_module
  use normal_mode_module
  use bands_module
  implicit none
end module

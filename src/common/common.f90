! ======================================================================
! Provides higher-level utilities, specific to the function of the code.
! Utilities include:
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module common_module
  use utils_module
  
  use structure_module
  use electronic_structure_module
  use bands_module
  use normal_mode_module
  use normal_mode_symmetry_module
  implicit none
end module

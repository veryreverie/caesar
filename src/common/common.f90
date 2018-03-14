! ======================================================================
! Provides higher-level utilities, specific to the function of the code.
! Utilities include:
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module common_module
  use utils_module
  
  use physical_constants_module
  use qpoints_module
  use atom_module
  use basic_symmetry_module
  use symmetry_module
  use structure_module
  use generate_qpoints_module
  use input_file_module
  use output_file_module
  use bands_module
  use normal_mode_module
  use normal_mode_symmetry_module
  implicit none
end module

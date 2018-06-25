! ======================================================================
! Interfaces to various electronic structure codes, and their files.
! ======================================================================
! This module is simply an interface for the various
!    electronic structure submodules.
module electronic_structure_module
  use electronic_structure_data_submodule
  use structure_file_submodule
  use quip_wrapper_submodule
  use electronic_structure_file_submodule
  use calculation_writer_submodule
  use calculation_runner_submodule
  use calculation_reader_submodule
  implicit none
end module

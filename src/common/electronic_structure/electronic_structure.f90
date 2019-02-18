! ======================================================================
! Interfaces to various electronic structure codes, and their files.
! ======================================================================
! This module is simply an interface for the various
!    electronic structure modules.
module electronic_structure_module
  use electronic_structure_data_module
  use quip_module
  
  use structure_file_module
  use electronic_structure_file_module
  use loto_splitting_module
  use kpoint_grid_module
  use calculation_writer_module
  use calculation_runner_module
  use calculation_reader_module
  use phonon_file_module
  
  use qe_wrapper_module, only: QeInputFile
  implicit none
end module

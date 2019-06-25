! ======================================================================
! Interfaces to various electronic structure codes, and their files.
! ======================================================================
! This module is simply an interface for the various
!    electronic structure modules.
module electronic_structure_module
  use electronic_structure_data_module
  use electronic_structure_common_module
  use quip_module
  use castep_module
  use qe_module
  use vasp_module
  
  use electronic_structure_file_module
  use calculation_writer_module
  use calculation_runner_module
  use calculation_reader_module
  implicit none
end module

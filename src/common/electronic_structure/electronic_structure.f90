! ======================================================================
! Interfaces to various electronic structure codes, and their files.
! ======================================================================
! This module is simply an interface for the various
!    electronic structure modules.
module caesar_electronic_structure_module
  use caesar_electronic_structure_data_module
  use caesar_electronic_structure_common_module
  use caesar_quip_module
  use caesar_castep_module
  use caesar_qe_module
  use caesar_vasp_module
  
  use caesar_electronic_structure_file_module
  use caesar_calculation_writer_module
  use caesar_calculation_runner_module
  use caesar_calculation_reader_module
  implicit none
end module

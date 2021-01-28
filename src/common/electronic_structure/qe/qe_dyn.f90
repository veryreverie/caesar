! ======================================================================
! Reads and writes QE dyn files.
! ======================================================================
module caesar_qe_dyn_module
  use caesar_utils_module
  use caesar_structure_module
  use caesar_normal_mode_module
  
  use caesar_dynamical_matrix_module
  implicit none
  
  private
  
  public :: read_qe_dynamical_matrix_file
contains

function read_qe_dynamical_matrix_file(directory,seedname,structure) &
   & result(output)
  ! TODO
end function
end module

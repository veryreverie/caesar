! ======================================================================
! Reads and writes QE dyn files.
! ======================================================================
module qe_dyn_module
  use utils_module
  use structure_module
  use normal_mode_module
  
  use dynamical_matrix_module
  implicit none
  
  private
  
  public :: read_qe_dynamical_matrix_file
contains

function read_qe_dynamical_matrix_file(directory,seedname,structure) &
   & result(output)
  ! TODO
end module

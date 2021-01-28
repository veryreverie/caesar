! ======================================================================
! Various routines involving force constant matrices,
!    in various co-ordinate systems.
! ======================================================================
! This module is simply an interface for the various dynamical matrices
!    modules.
! N.B. force constants in cartesian co-ordinates are Hessians,
!    and force constants in q-point co-ordinates are dynamical matrices.
module caesar_dynamical_matrices_module
  use caesar_unique_directions_module
  use caesar_min_images_module
  use caesar_dynamical_matrix_module
  use caesar_calculate_dynamical_matrices_module
  implicit none
end module

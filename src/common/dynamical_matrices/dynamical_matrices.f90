! ======================================================================
! Various routines involving force constant matrices,
!    in various co-ordinate systems.
! ======================================================================
! This module is simply an interface for the various dynamical matrices
!    modules.
! N.B. force constants in cartesian co-ordinates are Hessians,
!    and force constants in q-point co-ordinates are dynamical matrices.
module dynamical_matrices_module
  use unique_directions_module
  use min_images_module
  use cartesian_hessian_module
  use dynamical_matrix_module
  implicit none
end module

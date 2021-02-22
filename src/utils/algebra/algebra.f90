! ======================================================================
! Various mathematical functionality, including:
!    - Sparse vectors and matrices.
!    - Integer fractions.
!    - Fractional phases.
!    - Linear algebra.
!    - Diagonalisation and QR decomposition.
!    - Integer arrays.
!    - Groups.
! ======================================================================
! This module is simply an interface for the various algebra modules.
module caesar_algebra_module
  use caesar_mathematical_constants_module
  use caesar_linear_algebra_module
  use caesar_sparse_algebra_module
  
  use caesar_algebra_utils_module
  use caesar_fraction_module
  use caesar_phase_module
  use caesar_fraction_algebra_module
  use caesar_qr_decomposition_module
  use caesar_hermitian_eigenstuff_module
  use caesar_orthonormal_module
  use caesar_unitary_eigenstuff_module
  use caesar_integer_arrays_module
  use caesar_group_module
  use caesar_integer_complex_module
  use caesar_fraction_complex_module
  use caesar_lanczos_module
  use caesar_newton_raphson_module
  use caesar_tests_module
  implicit none
end module

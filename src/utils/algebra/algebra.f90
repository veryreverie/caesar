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
module algebra_module
  use mathematical_constants_module
  use linear_algebra_module
  use sparse_algebra_module
  
  use algebra_utils_module
  use fraction_module
  use phase_module
  use fraction_algebra_module
  use qr_decomposition_module
  use hermitian_eigenstuff_module
  use orthonormal_module
  use unitary_eigenstuff_module
  use integer_arrays_module
  use operator_group_module
  use integer_complex_module
  use fraction_complex_module
  use lanczos_module
  use newton_raphson_module
  use tests_module
  implicit none
end module

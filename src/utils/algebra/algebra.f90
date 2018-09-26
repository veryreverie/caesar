! ======================================================================
! Various mathematical functionality, including:
!    - Integer fractions.
!    - Fractional phases.
!    - Linear algebra.
!    - Diagonalisation and QR decomposition.
!    - Integer arrays.
!    - Groups.
! ======================================================================
! This module is simply an interface for the various algebra submodules.
module algebra_module
  use mathematical_constants_submodule
  use linear_algebra_submodule
  use algebra_utils_submodule
  use fraction_submodule
  use phase_submodule
  use fraction_algebra_submodule
  use qr_decomposition_submodule
  use hermitian_eigenstuff_submodule
  use orthonormal_submodule
  use unitary_eigenstuff_submodule
  use integer_arrays_submodule
  use group_submodule
  use integer_complex_submodule
  use fraction_complex_submodule
  use lanczos_submodule
  use tests_submodule
  implicit none
end module

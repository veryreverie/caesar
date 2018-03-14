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
  use eigenstuff_submodule
  use integer_arrays_submodule
  use group_submodule
  implicit none
end module

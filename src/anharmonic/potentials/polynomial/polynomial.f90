! ======================================================================
! Anharmonic functionality specific to the polynomial representation.
! ======================================================================
! This module is simply an interface for the various
!    polynomial anharmonic modules.
module caesar_polynomial_module
  use caesar_polynomial_stress_module
  use caesar_polynomial_potential_module
  implicit none
contains
subroutine startup_polynomial()
  implicit none
  
  call startup_polynomial_stress()
  call startup_polynomial_potential()
end subroutine
end module

! ======================================================================
! Anharmonic functionality specific to the polynomial representation.
! ======================================================================
! This module is simply an interface for the various
!    polynomial anharmonic modules.
module caesar_polynomial_module
  use caesar_polynomial_stress_module
  use caesar_polynomial_potential_module
  implicit none
  
  interface
    module subroutine startup_polynomial()
    end subroutine
  end interface
end module

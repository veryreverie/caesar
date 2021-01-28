! ======================================================================
! A unified interface to the different potential representations.
! ======================================================================
! This module is simply an interface for the various potentials modules.
module caesar_potentials_module
  use caesar_polynomial_module
  
  use caesar_potential_example_module
  implicit none
contains
subroutine startup_potentials()
  implicit none
  
  call startup_polynomial()
  call startup_potential_example()
end subroutine
end module

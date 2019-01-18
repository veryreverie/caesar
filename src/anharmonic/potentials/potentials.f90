! ======================================================================
! A unified interface to the different potential representations.
! ======================================================================
! This module is simply an interface for the various potentials modules.
module potentials_module
  use polynomial_module
  
  use potential_example_module
  implicit none
contains
subroutine startup_potentials()
  implicit none
  
  call startup_polynomial()
  call startup_potential_example()
end subroutine
end module

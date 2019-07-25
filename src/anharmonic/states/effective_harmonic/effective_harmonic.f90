! ======================================================================
! Harmonic bases of states.
! ======================================================================
! This module is simply an interface for the various harmonic states modules.
module effective_harmonic_module
  use harmonic_states_module
  use harmonic_basis_module
  implicit none
contains
subroutine startup_effective_harmonic()
  implicit none
  
  call startup_harmonic_states()
  call startup_harmonic_basis()
end subroutine
end module

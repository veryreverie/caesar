! ======================================================================
! States in various representations.
! ======================================================================
! This module is simply an interface for the various states modules.
module states_module
  use monomial_state_module
  use harmonic_state_module
  use polynomial_state_module
  use braket_module
  use full_subspace_basis_and_states_module
  use subspace_basis_pointer_module
  use subspace_states_pointer_module
  use state_conversion_module
  implicit none
contains
end module

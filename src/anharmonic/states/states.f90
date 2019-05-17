! ======================================================================
! States in various representations.
! ======================================================================
! This module is simply an interface for the various states modules.
module states_module
  use monomial_state_1d_module
  use monomial_state_2d_module
  use monomial_state_module
  use harmonic_state_module
  use polynomial_state_module
  use full_subspace_wavefunctions_module
  use split_qpoints_wavefunctions_module
  use full_subspace_basis_and_states_module
  use split_qpoints_basis_and_states_module
  use state_conversion_module
  implicit none
contains
subroutine startup_states()
  implicit none
  
  call startup_monomial_state()
  call startup_harmonic_state()
  call startup_polynomial_state()
  call startup_full_subspace_wavefunctions()
  call startup_split_qpoints_wavefunctions()
  call startup_full_subspace_basis_and_states()
  call startup_split_qpoints_basis_and_states()
end subroutine
end module

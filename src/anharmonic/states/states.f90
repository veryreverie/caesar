! ======================================================================
! States in various representations.
! ======================================================================
! This module is simply an interface for the various states modules.
module states_module
  use effective_harmonic_module
  
  use monomial_state_1d_module
  use monomial_state_2d_module
  use harmonic_state_1d_module
  use harmonic_state_2d_module
  use monomial_state_real_module
  use monomial_state_complex_module
  use harmonic_state_real_module
  use harmonic_state_complex_module
  use full_subspace_wavefunctions_module
  use split_qpoints_wavefunctions_module
  use wavevector_state_module
  use wavevector_states_module
  use full_subspace_basis_module
  use split_qpoints_basis_module
  implicit none
contains
subroutine startup_states()
  implicit none
  
  call startup_harmonic_states()
  call startup_monomial_state_real()
  call startup_monomial_state_complex()
  call startup_harmonic_state_real()
  call startup_harmonic_state_complex()
  call startup_wavevector_state()
  call startup_wavevector_states()
  call startup_full_subspace_wavefunctions()
  call startup_split_qpoints_wavefunctions()
  call startup_full_subspace_basis()
  call startup_split_qpoints_basis()
end subroutine
end module

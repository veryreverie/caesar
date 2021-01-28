! ======================================================================
! Anharmonic functionality which is independent of the choice of potential
!    representation.
! ======================================================================
! This module is simply an interface for the various
!    common anharmonic modules.
module caesar_anharmonic_common_module
  use caesar_integration_module
  
  use caesar_interpolated_supercell_module
  use caesar_max_displacement_module
  use caesar_subspace_coupling_module
  use caesar_subspace_monomial_module
  use caesar_degenerate_symmetry_module
  use caesar_anharmonic_data_module
  use caesar_subspace_wavefunctions_module
  use caesar_stress_prefactors_module
  use caesar_sparse_monomial_module
  use caesar_subspace_state_module
  use caesar_subspace_braket_module
  use caesar_basis_state_module
  use caesar_basis_states_module
  use caesar_pulay_module
  use caesar_abstract_classes_module
  use caesar_braket_module
  use caesar_stress_data_module
  use caesar_potential_data_module
  implicit none
contains
end module

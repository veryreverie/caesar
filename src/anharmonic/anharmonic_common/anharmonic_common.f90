! ======================================================================
! Anharmonic functionality which is independent of the choice of potential
!    representation.
! ======================================================================
! This module is simply an interface for the various
!    common anharmonic modules.
module anharmonic_common_module
  use integration_module
  
  use interpolated_supercell_module
  use max_displacement_module
  use subspace_coupling_module
  use subspace_monomial_module
  use degenerate_symmetry_module
  use anharmonic_data_module
  use subspace_wavefunctions_module
  use stress_prefactors_module
  use sparse_monomial_module
  use subspace_state_module
  use subspace_braket_module
  use basis_state_module
  use basis_states_module
  use pulay_module
  use abstract_classes_module
  use braket_module
  use stress_data_module
  use potential_data_module
  implicit none
contains
end module

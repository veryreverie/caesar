! ======================================================================
! Anharmonic functionality which is independent of the choice of potential
!    representation.
! ======================================================================
! This module is simply an interface for the various
!    common anharmonic submodules.
module anharmonic_common_module
  use degeneracy_module
  use coupled_subspaces_module
  use subspace_monomial_module
  use coupled_modes_module
  use mode_monomial_module
  use degenerate_symmetry_module
  use single_mode_state_module
  use subspace_product_state_module
  use potential_module
  use potential_pointer_module
  implicit none
contains
end module

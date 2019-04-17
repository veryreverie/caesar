! ======================================================================
! Anharmonic functionality which is independent of the choice of potential
!    representation.
! ======================================================================
! This module is simply an interface for the various
!    common anharmonic modules.
module anharmonic_common_module
  use subspace_coupling_module
  use subspace_monomial_module
  use degenerate_symmetry_module
  use anharmonic_data_module
  use energy_spectrum_module
  use subspace_wavefunctions_module
  use abstract_classes_module
  use braket_module
  implicit none
contains
end module

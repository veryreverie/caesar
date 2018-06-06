! ======================================================================
! The StructureData type, which holds all information about the structure
!    being modelled.
! Also contains all subsidiary types, including:
!    - Atoms.
!    - q-points.
!    - symmetries.
!    - R-vectors and G-vectors.
!    - Supercell data.
! ======================================================================
! This module is simply an interface for the various structure submodules.
module structure_module
  use physical_constants_submodule
  use basic_structure_submodule
  use atom_submodule
  use qpoint_submodule
  use generate_qpoints_submodule
  use calculate_symmetry_submodule
  use symmetry_submodule
  use structure_submodule
  use supercell_submodule
  use cartesian_displacement_submodule
  use cartesian_force_submodule
  implicit none
end module

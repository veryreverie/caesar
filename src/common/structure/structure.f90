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
! This module is simply an interface for the various structure modules.
module caesar_structure_module
  use caesar_atom_module
  
  use caesar_physical_constants_module
  use caesar_qpoint_module
  use caesar_generate_qpoints_module
  use caesar_symmetry_module
  use caesar_structure_data_module
  use caesar_supercell_module
  implicit none
end module

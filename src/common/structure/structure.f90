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
module structure_module
  use atom_module
  
  use physical_constants_module
  use basic_structure_module
  use qpoint_module
  use generate_qpoints_module
  use symmetry_module
  use structure_data_module
  use supercell_module
  implicit none
end module

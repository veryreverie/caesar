module calculate_symmetry_helper_module
contains

! ----------------------------------------------------------------------
! A helper function for calculate_symmetry.sh
!
! Takes the output of cellsym as run on a structure file
! Extracts the symmetries operations
! Appends these symmetry operations to the original structure file
! ----------------------------------------------------------------------
subroutine calculate_symmetry_helper(symmetry_filename,structure_filename)
  use string_module
  use structure_module
  implicit none
  
  ! filenames
  type(String), intent(in) :: symmetry_filename ! The output of cellsym
  type(String), intent(in) :: structure_filename ! A structure file
  
  ! file contents
  type(StructureData) :: structure
  
  ! Read structure data without symmetries.
  ! n.b. supercell is not used here, so a dummy is provided.
  structure = read_structure_file(structure_filename)
  
  ! Read symmetry data
  call read_symmetry_file(structure,symmetry_filename)
  
  ! Write both to file
  call write_structure_file(structure,structure_filename)
end subroutine
end module

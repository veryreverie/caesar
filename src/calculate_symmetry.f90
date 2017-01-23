module calculate_symmetry_helper_module
contains

! ----------------------------------------------------------------------
! A helper function for calculate_symmetry.sh
!
! Takes the output of cellsym as run on a structure file
! Extracts the symmetries operations
! Appends these symmetry operations to the original structure file
! ----------------------------------------------------------------------
subroutine calculate_symmetry_helper(filenames)
  use string_module
  use structure_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  ! filenames
  type(String) :: symmetry_filename ! The output of cellsym
  type(String) :: structure_filename ! A structure file
  
  ! file contents
  type(StructureData) :: structure
  
  symmetry_filename = filenames(1)
  structure_filename = filenames(2)
  
  structure = read_structure_file(structure_filename)
  
  call read_symmetry_file(structure,symmetry_filename)
  
  call write_structure_file(structure,structure_filename)
end subroutine
end module

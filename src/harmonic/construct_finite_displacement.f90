module construct_finite_displacement_module
  implicit none
contains

subroutine construct_finite_displacement(args)
  use string_module
  use structure_module
  implicit none
  
  type(String), intent(in) :: args(:)
  
  ! Input variables
  integer :: atom
  integer :: disp
  
  ! file names
  type(String) :: structure_filename
  type(String) :: positive_filename
  type(String) :: negative_filename
  
  type(StructureData) :: positive_structure
  type(StructureData) :: negative_structure

  ! Read in arguments
  atom = int(args(1))
  disp = int(args(2))
  structure_filename = args(3)
  positive_filename = args(4)
  negative_filename = args(5)
  
  ! Read in structure
  positive_structure = read_structure_file(structure_filename)
  negative_structure = read_structure_file(structure_filename)
  
  positive_structure%atoms(disp,atom) = positive_structure%atoms(disp,atom) &
                                    & + 0.01d0
  negative_structure%atoms(disp,atom) = negative_structure%atoms(disp,atom) &
                                    & - 0.02d0
  
  call write_structure_file(positive_structure,positive_filename)
  call write_structure_file(negative_structure,negative_filename)
end subroutine
end module

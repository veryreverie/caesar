! ======================================================================
! Provides helper functions for testing StructureData.
! ======================================================================
module structure_test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Writes out structure.dat in a manner compatible with old Caesar.
! ----------------------------------------------------------------------
subroutine write_old_structure_file(structure,filename)
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: filename
  
  type(RealMatrix), allocatable :: rotations(:)
  
  integer :: structure_file
  integer :: i
  
  rotations = calculate_cartesian_rotations(structure)
  
  structure_file = open_write_file(filename)
  call print_line(structure_file, 'Lattice')
  call print_line(structure_file, structure%lattice)
  call print_line(structure_file, 'Atoms')
  do i=1,structure%no_atoms
    call print_line(structure_file, structure%species(i) //' '//&
                                  & structure%mass(i)    //' '//&
                                  & structure%atoms(i))
  enddo
  call print_line(structure_file, 'Symmetry')
  do i=1,structure%no_symmetries
    call print_line(structure_file, rotations(i))
    call print_line(structure_file, structure%translations(i))
  enddo
  call print_line(structure_file, 'End')
  close(structure_file)
end subroutine
end module

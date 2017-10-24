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
  use ofile_module
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: filename
  
  type(RealMatrix), allocatable :: rotations(:)
  
  type(OFile) :: structure_file
  
  integer :: i
  
  rotations = structure%calculate_cartesian_rotations()
  
  structure_file = filename
  call structure_file%print_line('Lattice')
  call structure_file%print_line(structure%lattice)
  call structure_file%print_line('Atoms')
  do i=1,structure%no_atoms
    call structure_file%print_line( structure%species(i) //' '//&
                                  & structure%mass(i)    //' '//&
                                  & structure%atoms(i))
  enddo
  call structure_file%print_line('Symmetry')
  do i=1,size(structure%symmetries)
    call structure_file%print_line(rotations(i))
    call structure_file%print_line(structure%symmetries(i)%translation)
  enddo
  call structure_file%print_line('End')
end subroutine
end module

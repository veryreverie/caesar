! ======================================================================
! Provides helper functions for testing StructureData.
! ======================================================================
module caesar_structure_test_module
  use caesar_common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Writes out structure.dat in a manner compatible with old Caesar.
! ----------------------------------------------------------------------
subroutine write_old_structure_file(structure,filename)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: filename
  
  type(OFile) :: structure_file
  
  integer :: i
  
  structure_file = OFile(filename)
  call structure_file%print_line('Lattice')
  call structure_file%print_lines(structure%lattice)
  call structure_file%print_line('Atoms')
  do i=1,structure%no_atoms
    call structure_file%print_line( structure%atoms(i)%species() //' '// &
                                  & structure%atoms(i)%mass()    //' '// &
                                  & structure%atoms(i)%cartesian_position())
  enddo
  call structure_file%print_line('Symmetry')
  do i=1,size(structure%symmetries)
    call structure_file%print_lines(structure%symmetries(i)%cartesian_tensor)
    call structure_file%print_line(structure%symmetries(i)%translation)
  enddo
  call structure_file%print_line('End')
end subroutine
end module

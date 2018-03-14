! ======================================================================
! Provides helper functions for testing StructureData.
! ======================================================================
module structure_test_module
  use common_module
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
  
  structure_file = filename
  call structure_file%print_line('Lattice')
  call structure_file%print_line(structure%lattice)
  call structure_file%print_line('Atoms')
  do i=1,structure%no_atoms
    call structure_file%print_line( structure%atoms(i)%species() //' '// &
                                  & structure%atoms(i)%mass()    //' '// &
                                  & structure%atoms(i)%cartesian_position())
  enddo
  call structure_file%print_line('Symmetry')
  do i=1,size(structure%symmetries)
    call structure_file%print_line(structure%symmetries(i)%cartesian_rotation)
    call structure_file%print_line(structure%symmetries(i)%translation)
  enddo
  call structure_file%print_line('End')
end subroutine
end module

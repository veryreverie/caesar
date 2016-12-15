module construct_finite_displacement_module
  implicit none
contains

subroutine construct_finite_displacement(filenames)
  use constants, only : dp
  use file_io,   only : open_read_file, open_write_file
  implicit none
  
  character(32), intent(in) :: filenames(:)
  
  ! Input variables
  integer :: atom,disp,no_atoms
  real(dp) :: lattice(3,3)
  real(dp),allocatable :: atoms(:,:),mass(:)
  character(2),allocatable :: species(:)
  
  ! file units
  integer :: disp_file              ! disp.dat
  integer :: super_lattice_file     ! super_lattice.dat
  integer :: super_equilibrium_file ! super_equilibrium.dat
  integer :: positive_file          ! positive/structure.dat
  integer :: negative_file          ! negative/structure.dat
  
  integer :: i ! loop index

  ! Read in displacement
  disp_file = open_read_file(filenames(1))
  read(disp_file,*) atom, disp
  close(disp_file)

  ! Read in structure
  super_lattice_file = open_read_file(filenames(2))
  read(super_lattice_file,*) lattice(1,:)
  read(super_lattice_file,*) lattice(2,:)
  read(super_lattice_file,*) lattice(3,:)
  close(super_lattice_file)
  
  super_equilibrium_file = open_read_file(filenames(3))
  read(super_equilibrium_file,*) no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms),species(no_atoms))
  do i=1,no_atoms 
    read(super_equilibrium_file,*)species(i),mass(i),atoms(i,:)
  enddo
  close(super_equilibrium_file)

  ! Generate distorted structure
  positive_file = open_write_file(filenames(4))
  negative_file = open_write_file(filenames(5))
  write(positive_file,*)'Lattice'
  write(negative_file,*)'Lattice'
  do i=1,3
    write(positive_file,*)lattice(i,:)
    write(negative_file,*)lattice(i,:)
  enddo ! i
  write(positive_file,*)'Atoms'
  write(negative_file,*)'Atoms'
  do i=1,no_atoms
    if(i==atom)then
      atoms(i,disp)=atoms(i,disp)+0.01d0
    endif ! i==atom
    write(positive_file,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  do i=1,no_atoms
    if(i==atom)then
      atoms(i,disp)=atoms(i,disp)-0.02d0
    endif ! i==atom
    write(negative_file,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  write(positive_file,*)'Symmetry'
  write(negative_file,*)'Symmetry'
  write(positive_file,*)'End'
  write(negative_file,*)'End'
  close(positive_file)
  close(negative_file)

end subroutine
end module

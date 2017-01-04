module combine_forces_module
  implicit none
contains

subroutine combine_forces(filenames)
  use constants, only : dp, eV_per_A_to_au
  use file_io,   only : open_read_file, open_write_file
  implicit none
  
  character(100), intent(in) :: filenames(:)
  
  real(dp),allocatable :: positive(:),negative(:)
  integer,allocatable :: atom1(:),disp1(:),atom2(:),disp2(:)
  integer :: i,no_atoms
  
  ! file units
  integer :: super_equilibrium_file
  integer :: positive_file
  integer :: negative_file
  integer :: forces_file

  ! Read in number of atoms
  super_equilibrium_file = open_read_file(filenames(1))
  read(super_equilibrium_file,*) no_atoms
  close(super_equilibrium_file)

  ! Read in force constants
  allocate(positive(no_atoms*3),negative(no_atoms*3))
  allocate(atom1(no_atoms*3),disp1(no_atoms*3))
  allocate(atom2(no_atoms*3),disp2(no_atoms*3))
  
  positive_file = open_read_file(filenames(2))
  negative_file = open_read_file(filenames(3))
  do i=1,no_atoms*3
    read(positive_file,*)atom1(i),disp1(i),atom2(i),disp2(i),positive(i)
    read(negative_file,*)atom1(i),disp1(i),atom2(i),disp2(i),negative(i)
  enddo ! i
  close(positive_file)
  close(negative_file)

  ! Write out combined forces constants
  forces_file = open_write_file(filenames(4))
  do i=1,no_atoms*3
    write(forces_file,*) atom1(i), &
                       & disp1(i), &
                       & atom2(i), &
                       & disp2(i), &
                       & (positive(i)-negative(i))/(2.d0*0.01d0)*eV_per_A_to_au
  enddo
  close(forces_file)

end subroutine
end module

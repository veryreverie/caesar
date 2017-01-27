module combine_forces_module
  implicit none
contains

subroutine combine_forces(filenames)
  use constants, only : dp, eV_per_A_to_au
  use file_module,   only : open_read_file, open_write_file
  use string_module
  use structure_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  real(dp),allocatable :: positive(:),negative(:)
  integer,allocatable :: atom1(:),disp1(:),atom2(:),disp2(:)
  integer :: i
  
  type(StructureData) :: structure
  
  ! file units
  integer :: positive_file
  integer :: negative_file
  integer :: forces_file
  
  ! Read in number of atoms
  structure = read_structure_file(filenames(1))

  ! Read in force constants
  allocate(positive(structure%no_modes),negative(structure%no_modes))
  allocate(atom1(structure%no_modes),disp1(structure%no_modes))
  allocate(atom2(structure%no_modes),disp2(structure%no_modes))
  
  positive_file = open_read_file(filenames(2))
  negative_file = open_read_file(filenames(3))
  do i=1,structure%no_modes
    read(positive_file,*)atom1(i),disp1(i),atom2(i),disp2(i),positive(i)
    read(negative_file,*)atom1(i),disp1(i),atom2(i),disp2(i),negative(i)
  enddo ! i
  close(positive_file)
  close(negative_file)

  ! Write out combined forces constants
  forces_file = open_write_file(filenames(4))
  do i=1,structure%no_modes
    write(forces_file,*) atom1(i), &
                       & disp1(i), &
                       & atom2(i), &
                       & disp2(i), &
                       & (positive(i)-negative(i))/(2.d0*0.01d0)*eV_per_A_to_au
  enddo
  close(forces_file)

end subroutine
end module

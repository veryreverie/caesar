module equilibrium_frac_module
  implicit none
contains

subroutine equilibrium_frac(filenames)
  use constants,      only : dp
  use file_module,        only : open_read_file, open_write_file
  use linear_algebra, only : inv_33
  use string_module
  use structure_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  type(StructureData)   :: structure
  integer               :: i
  real(dp), allocatable :: frac_atoms(:,:)
  
  ! file names
  type(String) :: structure_filename
  type(String) :: super_equilibrium_frac_filename
  
  ! file units
  integer :: super_equilibrium_frac_file
  
  ! Read inputs
  structure_filename = filenames(1)
  super_equilibrium_frac_filename = filenames(2)
  
  ! Read in structure
  structure = read_structure_file(structure_filename)
  
  frac_atoms = matmul(structure%recip_lattice,structure%atoms)

  ! Write out fractional atomic positions
  super_equilibrium_frac_file = open_write_file(filenames(3))
  write(super_equilibrium_frac_file,*) structure%no_atoms
  do i=1,structure%no_atoms
    write(super_equilibrium_frac_file,*) structure%species(i), &
                                       & structure%mass(i),    &
                                       & frac_atoms(i,:)
  enddo
  close(super_equilibrium_frac_file)
end subroutine
end module

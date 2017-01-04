module equilibrium_frac_module
  implicit none
contains

subroutine equilibrium_frac(filenames)
  use constants,      only : dp
  use file_io,        only : open_read_file, open_write_file
  use linear_algebra, only : inv_33
  implicit none
  
  character(32), intent(in) :: filenames(:)
  
  integer :: i,no_atoms
  real(dp) :: lattice(3,3),trans_lattice(3,3),inv_lattice(3,3)
  real(dp),allocatable :: atoms(:,:),mass(:),frac_atoms(:,:)
  character(2),allocatable :: species(:)
  
  ! file units
  integer :: super_equilibrium_file
  integer :: super_lattice_file
  integer :: super_equilibrium_frac_file

  ! Read in atomic positions
  super_equilibrium_file = open_read_file(filenames(1))
  read(super_equilibrium_file,*)no_atoms
  allocate(atoms(no_atoms,3))
  allocate(mass(no_atoms))
  allocate(species(no_atoms))
  allocate(frac_atoms(no_atoms,3))
  do i=1,no_atoms
    read(super_equilibrium_file,*)species(i),mass(i),atoms(i,:)
  enddo
  close(super_equilibrium_file)

  ! Read in superlattice
  super_lattice_file = open_read_file(filenames(2))
  do i=1,3
    read(super_lattice_file,*)lattice(i,:)
  enddo
  close(super_lattice_file) 
  
  trans_lattice=transpose(lattice)
  call inv_33(trans_lattice,inv_lattice)
  do i=1,no_atoms
    frac_atoms(i,1:3) = atoms(i,1)*inv_lattice(1:3,1) &
                    & + atoms(i,2)*inv_lattice(1:3,2) &
                    & + atoms(i,3)*inv_lattice(1:3,3)
  enddo

  ! Write out fractional atomic positions
  super_equilibrium_frac_file = open_write_file(filenames(3))
  write(super_equilibrium_frac_file,*)no_atoms
  do i=1,no_atoms
    write(super_equilibrium_frac_file,*)species(i),mass(i),frac_atoms(i,:)
  enddo
  close(super_equilibrium_frac_file)
end subroutine
end module

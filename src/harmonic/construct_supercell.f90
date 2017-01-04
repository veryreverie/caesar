! Program to construct supercell from primitive cell

module construct_supercell_module
  implicit none
contains

subroutine construct_supercell(filenames)
  USE constants
  USE utils
  use linear_algebra,         only : determinant33
  use is_in_supercell_module, only : is_in_supercell
  use file_io,                only : open_read_file, open_write_file
  implicit none
  
  character(100), intent(in) :: filenames(:)
  
  ! Working variables
  INTEGER :: i,j,k,atom_counter,delta
  INTEGER :: dira,dirb,dirc
  REAL(dp) :: pos(3)
  ! Input variables
  INTEGER :: no_atoms,supercell(3,3),sc_size
  REAL(dp) :: lattice(3,3),super_lattice(3,3)
  REAL(dp),ALLOCATABLE :: atoms(:,:),mass(:),super_atoms(:,:),super_mass(:)
  CHARACTER(2),ALLOCATABLE :: species(:),super_species(:)
  
  ! file units
  integer :: lattice_file
  integer :: equilibrium_file
  integer :: supercell_file
  integer :: super_lattice_file
  integer :: super_equilibrium_file

  ! Read in lattice
  lattice_file = open_read_file(filenames(1))
  do i=1,3 
    read(lattice_file,*) lattice(i,1:3)
  enddo
  close(lattice_file)

  ! Read in equilibrium
  equilibrium_file = open_read_file(filenames(2))
  read(equilibrium_file,*) no_atoms
  allocate(atoms(no_atoms,3))
  allocate(mass(no_atoms))
  allocate(species(no_atoms))
  do i=1,no_atoms
    read(equilibrium_file,*) species(i),mass(i),atoms(i,1:3)
  enddo
  close(equilibrium_file)

  ! Read in supercell matrix
  supercell_file = open_read_file(filenames(3))
  do i=1,3
    read(supercell_file,*)supercell(i,1:3)
  enddo ! i
  close(supercell_file)
  
  sc_size=abs(determinant33(supercell))
  allocate(super_atoms(no_atoms*sc_size,3))
  allocate(super_mass(no_atoms*sc_size)) 
  allocate(super_species(no_atoms*sc_size))

  ! Generate supercell lattice
  super_lattice=MATMUL(supercell,lattice)
  super_lattice_file = open_write_file(filenames(4))
  do i=1,3
    write(super_lattice_file,*) super_lattice(i,1:3)
  enddo  
  close(super_lattice_file)

  ! Generate supercell atoms
  atom_counter=0
  delta=sc_size*3
  do i=1,no_atoms
    do dira=-delta,delta
      do dirb=-delta,delta
        do dirc=-delta,delta
          pos(:) = atoms(i,:)        &
               & + dira*lattice(1,:) &
               & + dirb*lattice(2,:) &
               & + dirc*lattice(3,:)
          if(is_in_supercell(pos,super_lattice))then
            atom_counter=atom_counter+1
            super_atoms(atom_counter,:)=pos(:)
            super_mass(atom_counter)=mass(i)
            super_species(atom_counter)=species(i)
          endif ! is in supercell          
        enddo ! dirc
      enddo ! dirb
    enddo ! dira
  enddo ! Loop over primitive cell atoms

  if(atom_counter/=(no_atoms*sc_size))then
    write(*,*)'Error in placing atoms in supercell! Please, try increasing delta.'
    STOP
  endif ! error in generating supercell 

  super_equilibrium_file = open_write_file(filenames(5))
  write(super_equilibrium_file,*) no_atoms*sc_size
  do i=1,no_atoms*sc_size
    write(super_equilibrium_file,*) super_species(i), &
                                  & super_mass(i),    &
                                  & super_atoms(i,:)
  enddo
  close(super_equilibrium_file)
end subroutine
end module

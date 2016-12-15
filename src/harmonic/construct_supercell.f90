! Program to construct supercell from primitive cell

module construct_supercell_module
  implicit none
contains

subroutine construct_supercell()
  USE constants
  USE utils
  use linear_algebra, only : determinant33
  use is_in_supercell_module, only : is_in_supercell
  IMPLICIT NONE
  ! Working variables
  INTEGER :: i,j,k,atom_counter,delta
  INTEGER :: dira,dirb,dirc
  REAL(dp) :: pos(3)
  ! Input variables
  INTEGER :: no_atoms,supercell(3,3),sc_size
  REAL(dp) :: lattice(3,3),super_lattice(3,3)
  REAL(dp),ALLOCATABLE :: atoms(:,:),mass(:),super_atoms(:,:),super_mass(:)
  CHARACTER(2),ALLOCATABLE :: species(:),super_species(:)
  

  ! Read in lattice
  open(1,file='lattice.dat')
  do i=1,3 
    read(1,*)lattice(i,1:3)
  enddo ! i
  close(1) ! lattice.dat

  ! Read in equilibrium
  open(1,file='equilibrium.dat')
  read(1,*)no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms))
  allocate(species(no_atoms))
  do i=1,no_atoms
    read(1,*)species(i),mass(i),atoms(i,1:3)
  enddo ! i
  close(1) ! equilibrium.dat


  ! Read in supercell matrix
  open(1,file='supercell.dat')
  do i=1,3
    read(1,*)supercell(i,1:3)
  enddo ! i
  close(1) ! supercell.dat
  sc_size=abs(determinant33(supercell))
  allocate(super_atoms(no_atoms*sc_size,3),super_mass(no_atoms*sc_size),&  
    &super_species(no_atoms*sc_size))

  ! Generate supercell lattice
  super_lattice=MATMUL(supercell,lattice)
  open(1,file='super_lattice.dat')
  do i=1,3
    write(1,*)super_lattice(i,1:3)
  enddo  
  close(1) ! super_lattice

  ! Generate supercell atoms
  atom_counter=0
  delta=sc_size*3
  do i=1,no_atoms
    do dira=-delta,delta
      do dirb=-delta,delta
        do dirc=-delta,delta
          pos(:)=atoms(i,:)+dira*lattice(1,:)+dirb*lattice(2,:)+dirc*lattice(3,:)
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

  open(1,file='super_equilibrium.dat')
  write(1,*)no_atoms*sc_size
  do i=1,no_atoms*sc_size
    write(1,*)super_species(i),super_mass(i),super_atoms(i,:)
  enddo ! i
  close(1)
end subroutine
end module

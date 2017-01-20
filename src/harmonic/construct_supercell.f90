! Program to construct supercell from primitive cell

module construct_supercell_module
  implicit none
contains

subroutine construct_supercell(filenames)
  use constants,        only : dp
  use linear_algebra,   only : determinant33, inv_33
  use file_io
  use structure_module
  use string_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  ! Parameters
  real(dp), parameter :: tol=1.d-2
  
  ! Working variables
  INTEGER :: i,atom_counter,delta
  INTEGER :: dira,dirb,dirc
  REAL(dp) :: pos(3)
  
  ! Input variables
  INTEGER :: supercell(3,3),sc_size
  REAL(dp) :: super_lattice(3,3)
  real(dp) :: inv_super_lattice(3,3)
  real(dp) :: frac_pos(3)
  REAL(dp),ALLOCATABLE :: super_atoms(:,:),super_mass(:)
  CHARACTER(2),ALLOCATABLE :: super_species(:)
  
  type(StructureData) :: structure
  type(StructureData) :: superstructure
  
  integer :: super_no_atoms
  
  ! filenames
  type(String) :: structure_filename
  type(String) :: supercell_filename
  type(String) :: superstructure_filename
  
  ! file units
  integer :: supercell_file
  
  ! Read filenames from input
  structure_filename = filenames(1)
  supercell_filename = filenames(2)
  superstructure_filename = filenames(3)
  
  ! Read in structure data
  structure = read_structure_file(structure_filename)
  
  ! Read in supercell matrix
  supercell_file = open_read_file(supercell_filename)
  do i=1,3
    read(supercell_file,*) supercell(i,1:3)
  enddo
  close(supercell_file)
  
  sc_size = abs(determinant33(supercell))
  super_no_atoms = structure%no_atoms*sc_size

  allocate(super_atoms(3,super_no_atoms))
  allocate(super_mass(super_no_atoms)) 
  allocate(super_species(super_no_atoms))

  ! Generate supercell lattice
  super_lattice = matmul(supercell,structure%lattice)

  ! Generate supercell atoms
  atom_counter=0
  delta=sc_size*3
  do i=1,structure%no_atoms
    do dira=-delta,delta
      do dirb=-delta,delta
        do dirc=-delta,delta
          pos(:) = structure%atoms(i,:)        &
               & + dira*structure%lattice(1,:) &
               & + dirb*structure%lattice(2,:) &
               & + dirc*structure%lattice(3,:)
          inv_super_lattice = inv_33(transpose(super_lattice))
          frac_pos = matmul(inv_super_lattice,pos)
          if ( frac_pos(1)>-tol .and. frac_pos(1)<(1.d0-tol) .and. &
             & frac_pos(2)>-tol .and. frac_pos(2)<(1.d0-tol) .and. &
             & frac_pos(3)>-tol .and. frac_pos(3)<(1.d0-tol) ) then
            atom_counter=atom_counter+1
            super_atoms(:,atom_counter)=pos(:)
            super_mass(atom_counter)=structure%mass(i)
            super_species(atom_counter)=structure%species(i)
          endif ! is in supercell          
        enddo ! dirc
      enddo ! dirb
    enddo ! dira
  enddo ! Loop over primitive cell atoms

  if(atom_counter/=(super_no_atoms))then
    write(*,*)'Error in placing atoms in supercell! Please, try increasing delta.'
    STOP
  endif ! error in generating supercell 
  
  ! Write output
  call new(superstructure,super_no_atoms,0)
  superstructure%lattice = super_lattice
  superstructure%species = super_species
  superstructure%mass = super_mass
  superstructure%atoms = super_atoms
  
  call write_structure_file(structure,superstructure_filename)
end subroutine
end module

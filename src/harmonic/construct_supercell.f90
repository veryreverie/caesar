! Program to construct supercell from primitive cell

module construct_supercell_module
  implicit none
contains

subroutine construct_supercell(filenames)
  use constants,        only : dp
  use linear_algebra,   only : determinant33, inv_33
  use file_module
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
  real(dp) :: frac_pos(3)
  
  type(StructureData) :: structure
  type(StructureData) :: structure_sc
  
  integer :: no_atoms_sc
  
  ! filenames
  type(String) :: structure_filename
  type(String) :: supercell_filename
  type(String) :: structure_sc_filename
  
  ! file units
  integer :: supercell_file
  
  ! Read filenames from input
  structure_filename = filenames(1)
  supercell_filename = filenames(2)
  structure_sc_filename = filenames(3)
  
  ! Read in structure data
  structure = read_structure_file(structure_filename)
  
  ! Read in supercell matrix
  supercell_file = open_read_file(supercell_filename)
  do i=1,3
    read(supercell_file,*) supercell(i,1:3)
  enddo
  close(supercell_file)
  
  sc_size = abs(determinant33(supercell))
  no_atoms_sc = structure%no_atoms*sc_size

  allocate(structure_sc%atoms(3,no_atoms_sc))
  allocate(structure_sc%mass(no_atoms_sc)) 
  allocate(structure_sc%species(no_atoms_sc))
  call new(structure_sc,no_atoms_sc,0)

  ! Generate supercell lattice
  structure_sc%lattice = matmul(supercell,structure%lattice)
  structure_sc%recip_lattice = inv_33(transpose(structure_sc%lattice))

  ! Generate supercell atoms
  atom_counter=0
  delta=sc_size*3
  do i=1,structure%no_atoms
    do dira=-delta,delta
      do dirb=-delta,delta
        do dirc=-delta,delta
          pos(:) = structure%atoms(i,:)        &
               & + matmul((/dira,dirb,dirc/),structure%lattice)
          frac_pos = matmul(structure%recip_lattice,pos)
          if (all(frac_pos>-tol) .and. all(frac_pos<1.d0-tol)) then
            atom_counter=atom_counter+1
            structure_sc%atoms(:,atom_counter)=pos(:)
            structure_sc%mass(atom_counter)=structure%mass(i)
            structure_sc%species(atom_counter)=structure%species(i)
          endif ! is in supercell          
        enddo ! dirc
      enddo ! dirb
    enddo ! dira
  enddo ! Loop over primitive cell atoms

  if(atom_counter/=(no_atoms_sc))then
    write(*,*)'Error in placing atoms in supercell! Please, try increasing delta.'
    STOP
  endif ! error in generating supercell 
  
  ! Write output
  call write_structure_file(structure,structure_sc_filename)
end subroutine
end module

! Program to construct supercell from primitive cell

module construct_supercell_module
  implicit none
contains

function construct_supercell(structure,supercell) result(structure_sc)
  use constants,        only : dp
  use linear_algebra,   only : determinant33, inv_33
  use file_module
  use structure_module
  use string_module
  implicit none
  
  ! Inputs
  type(StructureData), intent(in) :: structure
  integer,             intent(in) :: supercell(3,3)
  type(StructureData)             :: structure_sc
  
  ! Parameters
  real(dp), parameter :: tol=1.d-2
  
  ! Working variables
  INTEGER :: i,atom_counter,delta
  INTEGER :: dira,dirb,dirc
  REAL(dp) :: pos(3)
  
  ! Input variables
  INTEGER :: sc_size
  real(dp) :: frac_pos(3)
  
  integer :: no_atoms_sc
  
  sc_size = abs(determinant33(supercell))
  no_atoms_sc = structure%no_atoms*sc_size

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
  
end function
end module

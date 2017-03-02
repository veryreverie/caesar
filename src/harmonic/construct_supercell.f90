! Program to place atoms in the supercell.

module construct_supercell_module
  implicit none
contains

function construct_supercell(structure,supercell) result(structure_sc)
  use constants,        only : dp
  use linear_algebra,   only : determinant, invert, invert_int
  use file_module
  use structure_module
  use string_module
  use supercell_module
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(SupercellData), intent(in) :: supercell
  type(StructureData)             :: structure_sc
  
  ! Atomic positions.
  real(dp) :: atom_pos_sc(3)   ! Scaled fractional supercell co-ordinates.
  real(dp) :: copy_pos_cart(3) ! Cartesian co-ordinates.
  
  ! G-vector positions.
  integer :: gvector_sc(3)   ! Scaled fractional supercell co-ordinates.
  integer :: gvector_prim(3) ! Fractional primitive cell co-ordinates.
  
  ! Temporary variables
  integer :: i,j
  integer :: atom_counter
  integer :: no_atoms_sc
  
  no_atoms_sc = structure%no_atoms*supercell%sc_size

  call new(structure_sc,no_atoms_sc,0,supercell%sc_size)
  
  structure_sc%supercell = supercell

  ! Generate supercell lattice.
  structure_sc%lattice = matmul(supercell%supercell,structure%lattice)
  structure_sc%recip_lattice = invert(transpose(structure_sc%lattice))
  
  ! Generate atomic positions.
  do i=1,structure%no_atoms
    ! N.B. fractional supercell coordinates are scaled so that the supercell
    !    is a cartesian cube with side length sc_size.
    
    ! Calculate the position of atom i in scaled fractional supercell co-ords.
    atom_pos_sc = matmul(structure_sc%recip_lattice,structure%atoms(:,i)) &
                & * supercell%sc_size
    
    ! Loop accross G-vectors in the supercell, creating a copy of the atom i in
    !    the lattice copy corresponding to each G-vector.
    do j=1,supercell%sc_size
      
      gvector_prim = supercell%gvectors(:,j)
      
      ! Convert the G-vector from fractional primitive cell co-ordinates into
      !    scaled fractional supercell co-ordinates.
      gvector_sc = matmul(supercell%recip_supercell,gvector_prim)
      
      ! Translate the G-vector by supercell lattice vectors, s.t. the copy
      !    lies inside the primitive supercell.
      gvector_sc = modulo(gvector_sc+floor(atom_pos_sc),supercell%sc_size) &
               & - floor(atom_pos_sc)
      
      ! Convert the G-vector back into fractional primitive cell co-ordinates.
      gvector_prim = matmul(transpose(supercell%supercell),gvector_sc) &
                 & / supercell%sc_size
      
      ! Calculate the position of the copy in cartesian co-ordinates.
      copy_pos_cart = structure%atoms(:,i) &
                  & + matmul(transpose(structure%lattice),gvector_prim)
      
      ! Add the copy to the supercell.
      atom_counter = (j-1)*structure%no_atoms + i
      structure_sc%atoms(:,atom_counter) = copy_pos_cart
      structure_sc%mass(atom_counter) = structure%mass(i)
      structure_sc%species(atom_counter) = structure%species(i)
    enddo
  enddo
end function
end module

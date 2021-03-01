submodule (caesar_basic_structure_module) caesar_basic_structure_submodule
  use caesar_atom_module
contains

module procedure new_BasicStructure
  this%lattice_matrix = lattice_matrix
  this%atoms = atoms
end procedure

module procedure new_BasicStructure_constructed
  integer :: i,ialloc
  
  if (size(species)/=size(masses)) then
    call print_line(ERROR//': The number of atoms and atomic species do not &
       &match.')
    call err()
  elseif (size(masses)/=size(cartesian_positions)) then
    call print_line(ERROR//': The number of atomic species and positions do &
       &not match.')
    call err()
  endif
  
  this%lattice_matrix = lattice_matrix
  allocate(this%atoms(size(species)), stat=ialloc); call err(ialloc)
  do i=1,size(this%atoms)
    this%atoms(i) = BasicAtom(species(i), masses(i), cartesian_positions(i))
  enddo
end procedure

module procedure new_BasicSupercell
  if (abs(determinant(supercell_matrix))/=size(rvectors)) then
    call print_line(ERROR//': The size of the supercell, |S|, does not match &
       &the number of R-vectors provided.')
    call err()
  elseif (size(rvectors)/=size(gvectors)) then
    call print_line(ERROR//': The number of R-vectors and G-vectors do not &
       &match.')
    call err()
  elseif (size(atom_rvector_ids)/=size(atom_prim_ids)) then
    call print_line(ERROR//': The number of atom ids do not match.')
    call err()
  endif
  
  output%supercell_matrix = supercell_matrix
  output%rvectors         = rvectors
  output%gvectors         = gvectors
  output%atom_rvector_ids = atom_rvector_ids
  output%atom_prim_ids    = atom_prim_ids
end procedure

module procedure volume_BasicStructure
  output = abs(determinant(this%lattice_matrix))
end procedure
end submodule

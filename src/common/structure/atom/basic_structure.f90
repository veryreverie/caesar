! ======================================================================
! A minimal representation of the Structure class.
! ======================================================================
module caesar_basic_structure_module
  use caesar_utils_module
  
  use caesar_basic_atoms_module
  implicit none
  
  private
  
  public :: BasicStructure
  public :: BasicSupercell
  
  type, extends(NoDefaultConstructor) :: BasicStructure
    type(RealMatrix)             :: lattice_matrix
    type(BasicAtom), allocatable :: atoms(:)
  contains
    procedure, public :: volume => volume_BasicStructure
  end type
  
  interface BasicStructure
    module procedure new_BasicStructure
    module procedure new_BasicStructure_constructed
  end interface
  
  type, extends(NoDefaultConstructor) :: BasicSupercell
    type(IntMatrix)              :: supercell_matrix
    type(IntVector), allocatable :: rvectors(:)
    type(IntVector), allocatable :: gvectors(:)
    integer,         allocatable :: atom_rvector_ids(:)
    integer,         allocatable :: atom_prim_ids(:)
  end type
  
  interface BasicSupercell
    module procedure new_BasicSupercell
  end interface
contains
function new_BasicStructure(lattice_matrix,atoms) result(this)
  implicit none
  
  type(RealMatrix), intent(in) :: lattice_matrix
  type(BasicAtom),  intent(in) :: atoms(:)
  type(BasicStructure)         :: this
  
  this%lattice_matrix = lattice_matrix
  this%atoms = atoms
end function

function new_BasicStructure_constructed(lattice_matrix,species,masses, &
   & cartesian_positions) result(this)
  implicit none
  
  type(RealMatrix), intent(in) :: lattice_matrix
  type(String),     intent(in) :: species(:)
  real(dp),         intent(in) :: masses(:)
  type(RealVector), intent(in) :: cartesian_positions(:)
  type(BasicStructure)         :: this
  
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
end function

function new_BasicSupercell(supercell_matrix,rvectors,gvectors, &
   & atom_rvector_ids,atom_prim_ids) result(output)
  implicit none
  
  type(IntMatrix), intent(in) :: supercell_matrix
  type(IntVector), intent(in) :: rvectors(:)
  type(IntVector), intent(in) :: gvectors(:)
  integer,         intent(in) :: atom_rvector_ids(:)
  integer,         intent(in) :: atom_prim_ids(:)
  type(BasicSupercell)        :: output
  
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
end function

impure elemental function volume_BasicStructure(this) result(output)
  implicit none
  
  class(BasicStructure), intent(in) :: this
  real(dp)                          :: output
  
  output = abs(determinant(this%lattice_matrix))
end function
end module

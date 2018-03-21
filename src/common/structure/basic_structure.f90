! ======================================================================
! A minimal representation of the Structure class.
! ======================================================================
module basic_structure_submodule
  use utils_module
  implicit none
  
  private
  
  public :: BasicAtom
  public :: BasicStructure
  public :: BasicSymmetry
  public :: BasicSupercell
  
  type, extends(NoDefaultConstructor) :: BasicAtom
    type(String)     :: species
    real(dp)         :: mass
    type(RealVector) :: cartesian_position
  end type
  
  interface BasicAtom
    module procedure new_BasicAtom
  end interface
  
  type, extends(NoDefaultConstructor) :: BasicStructure
    type(RealMatrix)             :: lattice_matrix
    type(BasicAtom), allocatable :: atoms(:)
  end type
  
  interface BasicStructure
    module procedure new_BasicStructure
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
  
  type, extends(NoDefaultConstructor) :: BasicSymmetry
    type(IntMatrix)  :: rotation
    type(RealVector) :: translation
  end type
  
  interface BasicSymmetry
    module procedure new_BasicSymmetry
  end interface
contains

function new_BasicAtom(species,mass,cartesian_position) result(output)
  implicit none
  
  type(String),     intent(in) :: species
  real(dp),         intent(in) :: mass
  type(RealVector), intent(in) :: cartesian_position
  type(BasicAtom)              :: output
  
  output%species            = species
  output%mass               = mass
  output%cartesian_position = cartesian_position
end function

function new_BasicStructure(lattice_matrix,species,masses, &
   & cartesian_positions) result(output)
  implicit none
  
  type(RealMatrix), intent(in) :: lattice_matrix
  type(String),     intent(in) :: species(:)
  real(dp),         intent(in) :: masses(:)
  type(RealVector), intent(in) :: cartesian_positions(:)
  type(BasicStructure)         :: output
  
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
  
  output%lattice_matrix = lattice_matrix
  allocate(output%atoms(size(species)), stat=ialloc); call err(ialloc)
  do i=1,size(output%atoms)
    output%atoms(i) = BasicAtom(species(i), masses(i), cartesian_positions(i))
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

function new_BasicSymmetry(rotation,translation) result(output)
  implicit none
  
  type(IntMatrix),  intent(in) :: rotation
  type(RealVector), intent(in) :: translation
  
  type(BasicSymmetry) :: output
  
  output%rotation    = rotation
  output%translation = translation
end function
end module

! ======================================================================
! Information about an atom.
! ======================================================================
module caesar_atom_data_module
  use caesar_utils_module
  
  use caesar_basic_atoms_module
  implicit none
  
  private
  
  public :: AtomData
  public :: BasicAtom
  
  type :: AtomData
    ! Atom data.
    type(String),     private :: species_
    real(dp),         private :: mass_
    type(RealVector), private :: fractional_position_
    type(RealVector), private :: cartesian_position_
    
    ! Id numbers.
    integer, private :: id_
    integer, private :: prim_id_
    integer, private :: rvec_id_
    
    ! Copies of the structure's lattice and recip lattice.
    type(RealMatrix), private :: lattice_
    type(RealMatrix), private :: recip_lattice_
  contains
    ! Getters.
    procedure, public :: species
    procedure, public :: mass
    procedure, public :: fractional_position
    procedure, public :: cartesian_position
    
    procedure, public :: id
    procedure, public :: prim_id
    procedure, public :: rvec_id
    
    ! Setters.
    procedure, public :: set_species
    procedure, public :: set_mass
    procedure, public :: set_fractional_position
    procedure, public :: set_cartesian_position
  end type
  
  interface AtomData
    module procedure new_AtomData
  end interface
  
  interface BasicAtom
    module procedure new_BasicAtom_AtomData
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_AtomData(basic_atom,lattice,recip_lattice,id,prim_id,rvec_id) &
   & result(this)
  implicit none
  
  type(BasicAtom),  intent(in) :: basic_atom
  type(RealMatrix), intent(in) :: lattice
  type(RealMatrix), intent(in) :: recip_lattice
  integer,          intent(in) :: id
  integer,          intent(in) :: prim_id
  integer,          intent(in) :: rvec_id
  type(AtomData)               :: this
  
  this%lattice_ = lattice
  this%recip_lattice_ = recip_lattice
  
  this%species_ = basic_atom%species
  this%mass_ = basic_atom%mass
  
  call this%set_cartesian_position(basic_atom%cartesian_position)
  
  this%id_ = id
  this%prim_id_ = prim_id
  this%rvec_id_ = rvec_id
end function

! ----------------------------------------------------------------------
! Conversion to BasicAtom.
! ----------------------------------------------------------------------
impure elemental function new_BasicAtom_AtomData(this) result(output)
  implicit none
  
  type(AtomData), intent(in) :: this
  type(BasicAtom)            :: output
  
  output = BasicAtom(this%species(), this%mass(), this%cartesian_position())
end function

! ----------------------------------------------------------------------
! Getters.
! ----------------------------------------------------------------------
impure elemental function species(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  type(String)                :: output
  
  output = this%species_
end function

impure elemental function mass(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  real(dp)                    :: output
  
  output = this%mass_
end function

impure elemental function fractional_position(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  type(RealVector)            :: output
  
  output = this%fractional_position_
end function

impure elemental function cartesian_position(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  type(RealVector)            :: output
  
  output = this%cartesian_position_
end function

! The id of the atom in the supercell, i.e. the atom is supercell%atoms(id).
impure elemental function id(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  integer                     :: output
  
  output = this%id_
end function

! The id of the corresponding atom in the supercell's primitive cell.
impure elemental function prim_id(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  integer                     :: output
  
  output = this%prim_id_
end function

! The id of the R-vector which translates from the corresponding atom in the
!    supercell's primitive cell to this atom.
impure elemental function rvec_id(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  integer                     :: output
  
  output = this%rvec_id_
end function

! ----------------------------------------------------------------------
! Setters.
! ----------------------------------------------------------------------
subroutine set_species(this,species)
  implicit none
  
  class(AtomData), intent(inout) :: this
  type(String),    intent(in)    :: species
  
  this%species_ = species
end subroutine

subroutine set_mass(this,mass)
  implicit none
  
  class(AtomData), intent(inout) :: this
  real(dp),        intent(in)    :: mass
  
  this%mass_ = mass
end subroutine

subroutine set_fractional_position(this,fractional_position)
  implicit none
  
  class(AtomData),  intent(inout) :: this
  type(RealVector), intent(in)    :: fractional_position
  
  this%fractional_position_ = fractional_position
  this%cartesian_position_ = transpose(this%lattice_) * fractional_position
end subroutine

subroutine set_cartesian_position(this,cartesian_position)
  implicit none
  
  class(AtomData),  intent(inout) :: this
  type(RealVector), intent(in)    :: cartesian_position
  
  this%cartesian_position_ = cartesian_position
  this%fractional_position_ = this%recip_lattice_ * cartesian_position
end subroutine
end module

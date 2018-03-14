! ======================================================================
! Information about an atom.
! ======================================================================
module atom_module
  use utils_module
  implicit none
  
  private
  
  public :: AtomData
  
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
  
  type :: HasLattice
    type(RealMatrix) :: lattice
    type(RealMatrix) :: recip_lattice
  end type
  
  interface AtomData
    module procedure new_AtomData
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_AtomData(species,mass,cartesian_position,lattice,recip_lattice, &
   & id,prim_id,rvec_id) result(this)
  implicit none
  
  type(String),     intent(in) :: species
  real(dp),         intent(in) :: mass
  type(RealVector), intent(in) :: cartesian_position
  type(RealMatrix), intent(in) :: lattice
  type(RealMatrix), intent(in) :: recip_lattice
  integer,          intent(in) :: id
  integer,          intent(in) :: prim_id
  integer,          intent(in) :: rvec_id
  type(AtomData)               :: this
  
  this%lattice_ = lattice
  this%recip_lattice_ = recip_lattice
  
  this%species_ = species
  this%mass_ = mass
  
  call this%set_cartesian_position(cartesian_position)
  
  this%id_ = id
  this%prim_id_ = prim_id
  this%rvec_id_ = rvec_id
end function

! ----------------------------------------------------------------------
! Getters.
! ----------------------------------------------------------------------
function species(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  type(String)                :: output
  
  output = this%species_
end function

function mass(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  real(dp)                    :: output
  
  output = this%mass_
end function

function fractional_position(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  type(RealVector)            :: output
  
  output = this%fractional_position_
end function

function cartesian_position(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  type(RealVector)            :: output
  
  output = this%cartesian_position_
end function

! The id of the atom in the supercell, i.e. the atom is supercell%atoms(id).
function id(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  integer                     :: output
  
  output = this%id_
end function

! The id of the corresponding atom in the supercell's primitive cell.
function prim_id(this) result(output)
  implicit none
  
  class(AtomData), intent(in) :: this
  integer                     :: output
  
  output = this%prim_id_
end function

! The id of the R-vector which translates from the corresponding atom in the
!    supercell's primitive cell to this atom.
function rvec_id(this) result(output)
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

! ======================================================================
! Information about an atom.
! ======================================================================
module caesar_atom_data_module
  use caesar_utils_module
  
  use caesar_basic_atom_module
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
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    module function new_AtomData(basic_atom,lattice,recip_lattice,id, &
       & prim_id,rvec_id) result(this) 
      type(BasicAtom),  intent(in) :: basic_atom
      type(RealMatrix), intent(in) :: lattice
      type(RealMatrix), intent(in) :: recip_lattice
      integer,          intent(in) :: id
      integer,          intent(in) :: prim_id
      integer,          intent(in) :: rvec_id
      type(AtomData)               :: this
    end function
  end interface
  
  interface BasicAtom
    ! ----------------------------------------------------------------------
    ! Conversion to BasicAtom.
    ! ----------------------------------------------------------------------
    impure elemental module function new_BasicAtom_AtomData(this) &
       & result(output) 
      type(AtomData), intent(in) :: this
      type(BasicAtom)            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Getters.
    ! ----------------------------------------------------------------------
    impure elemental module function species(this) result(output) 
      class(AtomData), intent(in) :: this
      type(String)                :: output
    end function
  end interface
  
  interface
    impure elemental module function mass(this) result(output) 
      class(AtomData), intent(in) :: this
      real(dp)                    :: output
    end function
  end interface
  
  interface
    impure elemental module function fractional_position(this) result(output) 
      class(AtomData), intent(in) :: this
      type(RealVector)            :: output
    end function
  end interface
  
  interface
    impure elemental module function cartesian_position(this) result(output) 
      class(AtomData), intent(in) :: this
      type(RealVector)            :: output
    end function
  end interface
  
  interface
    ! The id of the atom in the supercell, i.e. the atom is supercell%atoms(id).
    impure elemental module function id(this) result(output) 
      class(AtomData), intent(in) :: this
      integer                     :: output
    end function
  end interface
  
  interface
    ! The id of the corresponding atom in the supercell's primitive cell.
    impure elemental module function prim_id(this) result(output) 
      class(AtomData), intent(in) :: this
      integer                     :: output
    end function
  end interface
  
  interface
    ! The id of the R-vector which translates from the corresponding atom in the
    !    supercell's primitive cell to this atom.
    impure elemental module function rvec_id(this) result(output) 
      class(AtomData), intent(in) :: this
      integer                     :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Setters.
    ! ----------------------------------------------------------------------
    module subroutine set_species(this,species) 
      class(AtomData), intent(inout) :: this
      type(String),    intent(in)    :: species
    end subroutine
  end interface
  
  interface
    module subroutine set_mass(this,mass) 
      class(AtomData), intent(inout) :: this
      real(dp),        intent(in)    :: mass
    end subroutine
  end interface
  
  interface
    module subroutine set_fractional_position(this,fractional_position) 
      class(AtomData),  intent(inout) :: this
      type(RealVector), intent(in)    :: fractional_position
    end subroutine
  end interface
  
  interface
    module subroutine set_cartesian_position(this,cartesian_position) 
      class(AtomData),  intent(inout) :: this
      type(RealVector), intent(in)    :: cartesian_position
    end subroutine
  end interface
end module

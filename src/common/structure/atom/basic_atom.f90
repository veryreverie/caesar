! ======================================================================
! A minimal representation of the Atom class.
! ======================================================================
module caesar_basic_atom_module
  use caesar_utils_module
  implicit none
  
  private
  
  public :: BasicAtom
  
  type, extends(NoDefaultConstructor) :: BasicAtom
    type(String)     :: species
    real(dp)         :: mass
    type(RealVector) :: cartesian_position
  end type
  
  interface BasicAtom
    module procedure new_BasicAtom
  end interface
  
  interface
    impure elemental module function new_BasicAtom(species,mass, &
       & cartesian_position) result(output) 
      type(String),     intent(in) :: species
      real(dp),         intent(in) :: mass
      type(RealVector), intent(in) :: cartesian_position
      type(BasicAtom)              :: output
    end function
  end interface
end module

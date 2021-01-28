! ======================================================================
! A minimal representation of the Atom class.
! ======================================================================
module caesar_basic_atoms_module
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
contains

impure elemental function new_BasicAtom(species,mass,cartesian_position) &
   & result(output)
  implicit none
  
  type(String),     intent(in) :: species
  real(dp),         intent(in) :: mass
  type(RealVector), intent(in) :: cartesian_position
  type(BasicAtom)              :: output
  
  output%species            = species
  output%mass               = mass
  output%cartesian_position = cartesian_position
end function
end module

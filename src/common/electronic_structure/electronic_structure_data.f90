! ======================================================================
! A class for holding calculated energy and forces.
! ======================================================================
module electronic_structure_data_submodule
  use utils_module
  
  use structure_module
  use normal_mode_module
  implicit none
  
  private
  
  public :: ElectronicStructure
  
  type, extends(Stringsable) :: ElectronicStructure
    real(dp)             :: energy
    type(CartesianForce) :: forces
  contains
    procedure, public :: read  => read_ElectronicStructure
    procedure, public :: write => write_ElectronicStructure
  end type
  
  interface ElectronicStructure
    module procedure new_ElectronicStructure
    module procedure new_ElectronicStructure_StringArray
  end interface
contains

! Constructor.
function new_ElectronicStructure(energy,forces) result(this)
  implicit none
  
  real(dp),             intent(in) :: energy
  type(CartesianForce), intent(in) :: forces
  type(ElectronicStructure)        :: this
  
  this%energy  = energy
  this%forces  = forces
end function

! I/O.
subroutine read_ElectronicStructure(this,input)
  implicit none
  
  class(ElectronicStructure), intent(out) :: this
  type(String),               intent(in)  :: input(:)
  
  real(dp)             :: energy
  type(CartesianForce) :: forces
  
  select type(this); type is(ElectronicStructure)
    energy = dble(input(2))
    forces = CartesianForce(StringArray(input(4:)))
    
    this = ElectronicStructure(energy,forces)
  end select
end subroutine

function write_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  select type(this); type is(ElectronicStructure)
    output = [ str('Energy (Hartree):'),      &
             & str(this%energy),              &
             & str('Forces (Hartree/Bohr):'), &
             & str(this%forces)               ]
  end select
end function

impure elemental function new_ElectronicStructure_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ElectronicStructure)     :: this
  
  this = input
end function
end module

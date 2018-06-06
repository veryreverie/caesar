! ======================================================================
! A class for holding calculated energy and forces.
! ======================================================================
module electronic_structure_data_submodule
  use utils_module
  
  use structure_module
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
  
  type(String), allocatable :: line(:)
  
  real(dp)             :: energy
  type(CartesianForce) :: forces
  
  line = split_line(input(1))
  energy = dble(line(3))
  forces = input(3:)
end subroutine

function write_ElectronicStructure(this) result(output)
  implicit none
  
  class(ElectronicStructure), intent(in) :: this
  type(String), allocatable              :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(2+size(this%forces)), stat=ialloc); call err(ialloc)
  output(1) = 'Energy (Hartree): '//this%energy
  output(2) = 'Forces (Hartree/Bohr):'
  output(3:) = str(this%forces)
end function

impure elemental function new_ElectronicStructure_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ElectronicStructure)     :: this
  
  this = input
end function
end module

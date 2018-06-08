! ======================================================================
! A force in complex mode co-ordinates.
! ======================================================================
module complex_mode_force_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_vector_submodule
  use complex_mode_vector_submodule
  implicit none
  
  private
  
  public :: ComplexModeForce
  
  type, extends(ComplexModeVector) :: ComplexModeForce
  contains
    procedure, public :: read  => read_ComplexModeForce
    procedure, public :: write => write_ComplexModeForce
  end type
  
  interface ComplexModeForce
    module procedure new_ComplexModeForce_ComplexModeVector
    module procedure new_ComplexModeForce_ComplexSingleModeVectors
    module procedure new_ComplexModeForce_StringArray
  end interface
contains

! Constructors.
function new_ComplexModeForce_ComplexModeVector(force) result(this)
  implicit none
  
  type(ComplexModeVector), intent(in) :: force
  type(ComplexModeForce)              :: this
  
  this%ComplexModeVector = force
end function

function new_ComplexModeForce_ComplexSingleModeVectors(forces) &
   & result(this)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: forces(:)
  type(ComplexModeForce)                    :: this
  
  this = ComplexModeForce(ComplexModeVector(forces))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeForce(this,input)
  implicit none
  
  class(ComplexModeForce), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeForce)
    this = ComplexModeForce(ComplexModeVector(StringArray(input)))
  end select
end subroutine

function write_ComplexModeForce(this) result(output)
  implicit none
  
  class(ComplexModeForce), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(ComplexModeForce)
    output = str(this%ComplexModeVector)
  end select
end function

impure elemental function new_ComplexModeForce_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexModeForce)        :: this
  
  this = input
end function
end module

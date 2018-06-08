! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module complex_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_vector_submodule
  use complex_mode_vector_submodule
  implicit none
  
  private
  
  public :: ComplexModeDisplacement
  
  type, extends(ComplexModeVector) :: ComplexModeDisplacement
  contains
    procedure, public :: read  => read_ComplexModeDisplacement
    procedure, public :: write => write_ComplexModeDisplacement
  end type
  
  interface ComplexModeDisplacement
    module procedure new_ComplexModeDisplacement_ComplexModeVector
    module procedure new_ComplexModeDisplacement_ComplexSingleModeVectors
    module procedure new_ComplexModeDisplacement_StringArray
  end interface
contains

! Constructors.
function new_ComplexModeDisplacement_ComplexModeVector(displacement) &
   & result(this)
  implicit none
  
  type(ComplexModeVector), intent(in) :: displacement
  type(ComplexModeDisplacement)       :: this
  
  this%ComplexModeVector = displacement
end function

function new_ComplexModeDisplacement_ComplexSingleModeVectors(displacements) &
   & result(this)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: displacements(:)
  type(ComplexModeDisplacement)             :: this
  
  this = ComplexModeDisplacement(ComplexModeVector(displacements))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeDisplacement(this,input)
  implicit none
  
  class(ComplexModeDisplacement), intent(out) :: this
  type(String),                   intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeDisplacement)
    this = ComplexModeDisplacement(ComplexModeVector(StringArray(input)))
  end select
end subroutine

function write_ComplexModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  select type(this); type is(ComplexModeDisplacement)
    output = str(this%ComplexModeVector)
  end select
end function

impure elemental function new_ComplexModeDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexModeDisplacement) :: this
  
  this = input
end function
end module

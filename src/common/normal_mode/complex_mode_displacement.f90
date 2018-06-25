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
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(ComplexModeVector) :: ComplexModeDisplacement
  contains
    procedure, public :: read  => read_ComplexModeDisplacement
    procedure, public :: write => write_ComplexModeDisplacement
  end type
  
  interface ComplexModeDisplacement
    module procedure new_ComplexModeDisplacement_ComplexModeVector
    module procedure new_ComplexModeDisplacement_ComplexSingleModeVectors
    module procedure new_ComplexModeDisplacement_ComplexModes
    module procedure new_ComplexModeDisplacement_StringArray
  end interface
  
  interface operator(*)
    module procedure multiply_real_ComplexModeDisplacement
    module procedure multiply_ComplexModeDisplacement_real
    module procedure multiply_complex_ComplexModeDisplacement
    module procedure multiply_ComplexModeDisplacement_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexModeDisplacement_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexModeDisplacement_ComplexModeDisplacement
  end interface
  
  interface sum
    module procedure sum_ComplexModeDisplacements
  end interface
  
  interface operator(-)
    module procedure negative_ComplexModeDisplacement
    module procedure subtract_ComplexModeDisplacement_ComplexModeDisplacement
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

function new_ComplexModeDisplacement_ComplexModes(modes,displacements) &
   & result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: modes(:)
  complex(dp),       intent(in) :: displacements(:)
  type(ComplexModeDisplacement) :: this
  
  this = ComplexModeDisplacement(ComplexModeVector(modes,displacements))
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_ComplexModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                      intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this*that%ComplexModeVector)
end function

impure elemental function multiply_ComplexModeDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  real(dp),                      intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this%ComplexModeVector*that)
end function

impure elemental function multiply_complex_ComplexModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  complex(dp),                   intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this*that%ComplexModeVector)
end function

impure elemental function multiply_ComplexModeDisplacement_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this%ComplexModeVector*that)
end function

impure elemental function divide_ComplexModeDisplacement_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(this%ComplexModeVector/that)
end function

impure elemental function add_ComplexModeDisplacement_ComplexModeDisplacement(&
   & this,that) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement( this%ComplexModeVector &
                                & + that%ComplexModeVector)
end function

function sum_ComplexModeDisplacements(this) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this(:)
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(sum(this%ComplexModeVector))
end function

impure elemental function negative_ComplexModeDisplacement(this) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement(-this%ComplexModeVector)
end function

impure elemental function                                                &
   & subtract_ComplexModeDisplacement_ComplexModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: that
  type(ComplexModeDisplacement)             :: output
  
  output = ComplexModeDisplacement( this%ComplexModeVector &
                                & - that%ComplexModeVector)
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

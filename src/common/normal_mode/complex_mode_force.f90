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
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
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
  
  interface operator(*)
    module procedure multiply_real_ComplexModeForce
    module procedure multiply_ComplexModeForce_real
    module procedure multiply_complex_ComplexModeForce
    module procedure multiply_ComplexModeForce_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexModeForce_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexModeForce_ComplexModeForce
  end interface
  
  interface sum
    module procedure sum_ComplexModeForces
  end interface
  
  interface operator(-)
    module procedure negative_ComplexModeForce
    module procedure subtract_ComplexModeForce_ComplexModeForce
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
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_ComplexModeForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),               intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this*that%ComplexModeVector)
end function

impure elemental function multiply_ComplexModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  real(dp),               intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%ComplexModeVector*that)
end function

impure elemental function multiply_complex_ComplexModeForce(this,that) &
   & result(output)
  implicit none
  
  complex(dp),            intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this*that%ComplexModeVector)
end function

impure elemental function multiply_ComplexModeForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  complex(dp),            intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%ComplexModeVector*that)
end function

impure elemental function divide_ComplexModeForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  complex(dp),            intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%ComplexModeVector/that)
end function

impure elemental function add_ComplexModeForce_ComplexModeForce(this, &
   & that) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%ComplexModeVector + that%ComplexModeVector)
end function

function sum_ComplexModeForces(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this(:)
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(sum(this%ComplexModeVector))
end function

impure elemental function negative_ComplexModeForce(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(-this%ComplexModeVector)
end function

impure elemental function subtract_ComplexModeForce_ComplexModeForce( &
   & this,that) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%ComplexModeVector - that%ComplexModeVector)
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

! ======================================================================
! A force in complex mode co-ordinates.
! ======================================================================
module complex_mode_force_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_force_submodule
  implicit none
  
  private
  
  public :: ComplexModeForce
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: ComplexModeForce
    type(ComplexSingleForce), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_ComplexModeForce
    procedure, public :: write => write_ComplexModeForce
  end type
  
  interface ComplexModeForce
    module procedure new_ComplexModeForce
    module procedure new_ComplexModeForce_ComplexModes
    module procedure new_ComplexModeForce_StringArray
  end interface
  
  interface size
    module procedure size_ComplexModeForce
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

! Constructors and size() function.
function new_ComplexModeForce(forces) result(this)
  implicit none
  
  type(ComplexSingleForce), intent(in) :: forces(:)
  type(ComplexModeForce)               :: this
  
  this%vectors = forces
end function

function new_ComplexModeForce_ComplexModes(modes,forces) result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: modes(:)
  complex(dp),       intent(in) :: forces(:)
  type(ComplexModeForce)        :: this
  
  this = ComplexModeForce(ComplexSingleForce(modes,forces))
end function

function size_ComplexModeForce(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  integer                            :: output
  
  output = size(this%vectors)
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
  
  output = ComplexModeForce(this*that%vectors)
end function

impure elemental function multiply_ComplexModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  real(dp),               intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%vectors*that)
end function

impure elemental function multiply_complex_ComplexModeForce(this,that) &
   & result(output)
  implicit none
  
  complex(dp),            intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this*that%vectors)
end function

impure elemental function multiply_ComplexModeForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  complex(dp),            intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%vectors*that)
end function

impure elemental function divide_ComplexModeForce_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  complex(dp),            intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(this%vectors/that)
end function

impure elemental function add_ComplexModeForce_ComplexModeForce(this, &
   & that) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  integer :: i,j
  
  output = this
  do i=1,size(that)
    j = first(this%vectors%id==that%vectors(i)%id, default=0)
    if (j==0) then
      output%vectors = [output%vectors, that%vectors(i)]
    else
      output%vectors(j) = output%vectors(j) + that%vectors(i)
    endif
  enddo
end function

function sum_ComplexModeForces(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this(:)
  type(ComplexModeForce)             :: output
  
  integer :: i
  
  if (size(this)==0) then
    call print_line(ERROR//': Trying to sum an empty list.')
    call err()
  endif
  
  output = this(1)
  do i=2,size(this)
    output = output + this(i)
  enddo
end function

impure elemental function negative_ComplexModeForce(this) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce)             :: output
  
  output = ComplexModeForce(-this%vectors)
end function

impure elemental function subtract_ComplexModeForce_ComplexModeForce( &
   & this,that) result(output)
  implicit none
  
  type(ComplexModeForce), intent(in) :: this
  type(ComplexModeForce), intent(in) :: that
  type(ComplexModeForce)             :: output
  
  output = this + (-that)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeForce(this,input)
  implicit none
  
  class(ComplexModeForce), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeForce)
    this = ComplexModeForce(ComplexSingleForce(input))
  end select
end subroutine

function write_ComplexModeForce(this) result(output)
  implicit none
  
  class(ComplexModeForce), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  select type(this); type is(ComplexModeForce)
    output = str(this%vectors)
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

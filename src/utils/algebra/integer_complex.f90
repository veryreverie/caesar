! ======================================================================
! Complex numbers with integer real and imaginary parts.
! ======================================================================
module integer_complex_submodule
  use io_module
  implicit none
  
  private
  
  public :: IntComplex
  public :: real
  public :: aimag
  public :: conjg
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  
  type, extends(Stringable) :: IntComplex
    integer, private :: real_
    integer, private :: imag_
  contains
    procedure, public :: to_String => to_String_IntComplex
  end type
  
  interface IntComplex
    module procedure new_IntComplex
  end interface
  
  interface real
    module procedure real_IntComplex
  end interface
  
  interface aimag
    module procedure aimag_IntComplex
  end interface
  
  interface conjg
    module procedure conjg_IntComplex
  end interface
  
  interface operator(+)
    module procedure add_IntComplex_IntComplex
    module procedure add_IntComplex_integer
    module procedure add_integer_IntComplex
  end interface
  
  interface operator(-)
    module procedure subtract_IntComplex_IntComplex
    module procedure subtract_IntComplex_integer
    module procedure subtract_integer_IntComplex
  end interface
  
  interface operator(*)
    module procedure multiply_IntComplex_IntComplex
    module procedure multiply_IntComplex_integer
    module procedure multiply_integer_IntComplex
  end interface
contains

! Constructor.
function new_IntComplex(real_part,imaginary_part) result(this)
  implicit none
  
  integer, intent(in)           :: real_part
  integer, intent(in), optional :: imaginary_part
  type(IntComplex)              :: this
  
  this%real_ = real_part
  if (present(imaginary_part)) then
    this%imag_ = imaginary_part
  else
    this%imag_ = 0
  endif
end function

! Complex functions (real, aimag and conjg).
function real_IntComplex(this) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  integer                      :: output
  
  output = this%real_
end function

function aimag_IntComplex(this) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  integer                      :: output
  
  output = this%imag_
end function

function conjg_IntComplex(this) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  type(IntComplex)             :: output
  
  output = IntComplex(real(this),-aimag(this))
end function

! Arithmetic.
impure elemental function add_IntComplex_IntComplex(this,that) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  type(IntComplex), intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end function

impure elemental function add_IntComplex_integer(this,that) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  integer,          intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = real(this) + that
  output%imag_ = aimag(this)
end function

impure elemental function add_integer_IntComplex(this,that) result(output)
  implicit none
  
  integer,          intent(in) :: this
  type(IntComplex), intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = this + real(that)
  output%imag_ =        aimag(that)
end function

impure elemental function subtract_IntComplex_IntComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  type(IntComplex), intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end function

impure elemental function subtract_IntComplex_integer(this,that) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  integer,          intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = real(this) - that
  output%imag_ = aimag(this)
end function

impure elemental function subtract_integer_IntComplex(this,that) result(output)
  implicit none
  
  integer,          intent(in) :: this
  type(IntComplex), intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = this - real(that)
  output%imag_ =      - aimag(that)
end function

impure elemental function multiply_IntComplex_IntComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  type(IntComplex), intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end function

impure elemental function multiply_IntComplex_integer(this,that) result(output)
  implicit none
  
  type(IntComplex), intent(in) :: this
  integer,          intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = real(this) * that
  output%imag_ = aimag(this) * that
end function

impure elemental function multiply_integer_IntComplex(this,that) result(output)
  implicit none
  
  integer,          intent(in) :: this
  type(IntComplex), intent(in) :: that
  type(IntComplex)             :: output
  
  output%real_ = this * real(that)
  output%imag_ = this * aimag(that)
end function

! I/O.
function to_String_IntComplex(this) result(output)
  implicit none
  
  class(IntComplex), intent(in) :: this
  type(String)                  :: output
  
  if (aimag(this)>=0) then
    output = real(this)//'+'//aimag(this)//'i'
  else
    output = real(this)//'-'//-aimag(this)//'i'
  endif
end function
end module

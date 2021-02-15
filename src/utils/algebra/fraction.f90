! ======================================================================
! Integer fractions in an exact representation.
! ======================================================================
module caesar_fraction_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_algebra_utils_module
  implicit none
  
  private
  
  public :: IntFraction
  public :: operator(==)
  public :: operator(/=)
  public :: operator(<)
  public :: operator(>)
  public :: operator(<=)
  public :: operator(>=)
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: int
  public :: dble
  public :: frac
  public :: is_int
  public :: modulo
  public :: abs
  
  type, extends(Stringable) :: IntFraction
    integer, private :: n_
    integer, private :: d_
  contains
    ! Getters
    procedure, public :: denominator
    procedure, public :: numerator
    
    ! Fraction simplification.
    procedure, private :: simplify
    
    ! I/O.
    procedure, public :: read  => read_IntFraction
    procedure, public :: write => write_IntFraction
  end type
  
  ! Comparison.
  interface operator(==)
    module procedure equality_IntFraction_IntFraction
    module procedure equality_IntFraction_integer
    module procedure equality_integer_IntFraction
  end interface
  
  interface operator(/=)
    module procedure non_equality_IntFraction_IntFraction
    module procedure non_equality_IntFraction_integer
    module procedure non_equality_integer_IntFraction
  end interface
  
  interface operator(<)
    module procedure lt_IntFraction_IntFraction
    module procedure lt_IntFraction_integer
    module procedure lt_integer_IntFraction
  end interface
  
  interface operator(>)
    module procedure gt_IntFraction_IntFraction
    module procedure gt_IntFraction_integer
    module procedure gt_integer_IntFraction
  end interface
  
  interface operator(<=)
    module procedure le_IntFraction_IntFraction
    module procedure le_IntFraction_integer
    module procedure le_integer_IntFraction
  end interface
  
  interface operator(>=)
    module procedure ge_IntFraction_IntFraction
    module procedure ge_IntFraction_integer
    module procedure ge_integer_IntFraction
  end interface
  
  ! Arithmetic.
  interface operator(+)
    module procedure add_IntFraction_IntFraction
    module procedure add_IntFraction_integer
    module procedure add_integer_IntFraction
  end interface
  
  interface operator(-)
    module procedure subtract_IntFraction_IntFraction
    module procedure subtract_IntFraction_integer
    module procedure subtract_integer_IntFraction
  end interface
  
  interface operator(*)
    module procedure multiply_IntFraction_IntFraction
    module procedure multiply_IntFraction_integer
    module procedure multiply_integer_IntFraction
  end interface
  
  interface operator(/)
    module procedure divide_IntFraction_IntFraction
    module procedure divide_IntFraction_integer
    module procedure divide_integer_IntFraction
  end interface
  
  ! Constructor.
  interface IntFraction
    module procedure new_IntFraction
    module procedure new_IntFraction_String
  end interface
  
  ! Conversions to and from other types.
  interface int
    module procedure int_IntFraction
  end interface
  
  interface dble
    module procedure dble_IntFraction
  end interface
  
  interface frac
    module procedure frac_character
    module procedure frac_String
    module procedure frac_integer
    module procedure frac_integers
  end interface
  
  ! Check whether or not the fraction is an integer.
  interface is_int
    module procedure is_int_IntFraction
  end interface
  
  ! The fraction modulo an integer.
  interface modulo
    module procedure modulo_IntFraction_integer
  end interface
  
  ! Negative.
  interface operator(-)
    module procedure negative_IntFraction
  end interface
  
  ! Absolute value.
  interface abs
    module procedure abs_IntFraction
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_IntFraction(numerator,denominator) result(this)
  implicit none
  
  integer, intent(in) :: numerator
  integer, intent(in) :: denominator
  type(IntFraction)   :: this
  
  this%n_ = numerator
  this%d_ = denominator
  call this%simplify()
end function

! ----------------------------------------------------------------------
! Getters.
! ----------------------------------------------------------------------
impure elemental function numerator(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer                        :: output
  
  output = this%n_
end function

impure elemental function denominator(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer                        :: output
  
  output = this%d_
end function

! ----------------------------------------------------------------------
! Simplifies the fraction.
! ----------------------------------------------------------------------
subroutine simplify(this)
  implicit none
  
  class(IntFraction), intent(inout) :: this
  
  integer :: common_denominator
  
  common_denominator = gcd(this%n_, this%d_)
  this%n_ = this%n_ / common_denominator
  this%d_ = this%d_ / common_denominator
  
  if (this%d_<0) then
    this%n_ = - this%n_
    this%d_ = - this%d_
  endif
end subroutine

! ----------------------------------------------------------------------
! Conversions to and from other types.
! ----------------------------------------------------------------------
! Conversion to integer.
! As with int(real), rounds down non-integer fractions.
impure elemental function int_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  integer                       :: output
  
  output = this%n_/this%d_
end function

! Conversion to real(dp).
impure elemental function dble_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  real(dp)                      :: output
  
  output = real(this%n_,dp) / this%d_
end function

! Conversion from character(*).
impure elemental function frac_character(input) result(output)
  implicit none
  
  character(*), intent(in) :: input
  type(IntFraction)        :: output
  
  output = frac(str(input))
end function

! Conversion from String.
impure elemental function frac_String(input) result(output)
  implicit none
  
  type(String), intent(in) :: input
  type(IntFraction)        :: output
  
  output = IntFraction(input)
end function

! Conversion from integer.
impure elemental function frac_integer(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  type(IntFraction)   :: output
  
  output = IntFraction(input,1)
end function

! Conversion from numerator and denominator.
function frac_integers(numerator,denominator) result(output)
  implicit none
  
  integer, intent(in) :: numerator
  integer, intent(in) :: denominator
  type(IntFraction)   :: output
  
  output = IntFraction(numerator,denominator)
end function

! ----------------------------------------------------------------------
! Comparison.
! ----------------------------------------------------------------------
impure elemental function equality_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this%n_==that%n_ .and. this%d_==that%d_
end function

impure elemental function equality_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = this%n_==that .and. this%d_==1
end function

impure elemental function equality_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this==that%n_ .and. that%d_==1
end function

impure elemental function non_equality_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

impure elemental function lt_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this%n_*that%d_ < that%n_*this%d_
end function

impure elemental function lt_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = this%n_ < that*this%d_
end function

impure elemental function lt_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this*that%d_ < that%n_
end function

impure elemental function gt_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this%n_*that%d_ > that%n_*this%d_
end function

impure elemental function gt_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = this%n_ > that*this%d_
end function

impure elemental function gt_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this*that%d_ > that%n_
end function

impure elemental function le_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this > that
end function

impure elemental function le_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = .not. this > that
end function

impure elemental function le_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this > that
end function

impure elemental function ge_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this < that
end function

impure elemental function ge_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = .not. this < that
end function

impure elemental function ge_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this < that
end function

! ----------------------------------------------------------------------
! Addition.
! ----------------------------------------------------------------------
! a/b + c/d = (ad+bc)/(bd).
impure elemental function add_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_+that%n_*this%d_, this%d_*that%d_)
end function

! a/b + c = (a+bc)/b.
impure elemental function add_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_+that*this%d_, this%d_)
end function

! a + b/c = (ac+b)/c.
impure elemental function add_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%d_+that%n_, that%d_)
end function

! ----------------------------------------------------------------------
! Subtraction.
! ----------------------------------------------------------------------
! a/b + c/d = (ad-bc)/(bd).
impure elemental function subtract_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_-that%n_*this%d_, this%d_*that%d_)
end function

! a/b - c = (a-bc)/b.
impure elemental function subtract_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_-that*this%d_, this%d_)
end function

! a - b/c = (ac-b)/c.
impure elemental function subtract_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%d_-that%n_, that%d_)
end function

! ----------------------------------------------------------------------
! Multiplication.
! ----------------------------------------------------------------------
! a/b * c/d = (ac)/(bd).
impure elemental function multiply_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%n_, this%d_*that%d_)
end function

! a/b * c = ac/b.
impure elemental function multiply_IntFraction_integer(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that, this%d_)
end function

! a * b/c = ab/c.
impure elemental function multiply_integer_IntFraction(this,that) &
   & result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%n_, that%d_)
end function

! ----------------------------------------------------------------------
! Division.
! ----------------------------------------------------------------------
! a/b / c/d = (ad)/(bc).
impure elemental function divide_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_, this%d_*that%n_)
end function

! a/b / c = a/(bc).
impure elemental function divide_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_, this%d_*that)
end function

! a / b/c = ac/b.
impure elemental function divide_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%d_, that%n_)
end function

! ----------------------------------------------------------------------
! Whether or not the IntFraction is an integer.
! ----------------------------------------------------------------------
! Equivalent to whether or not the denominator = 1.
impure elemental function is_int_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  logical                       :: output
  
  output = this%d_==1
end function

! ----------------------------------------------------------------------
! A fraction modulo an integer.
! ----------------------------------------------------------------------
impure elemental function modulo_IntFraction_integer(this,that) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  integer,           intent(in) :: that
  type(IntFraction)             :: output
  
  output = IntFraction(modulo(this%n_,this%d_*that), this%d_)
end function

! ----------------------------------------------------------------------
! Negative.
! ----------------------------------------------------------------------
impure elemental function negative_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  type(IntFraction)             :: output
  
  output = IntFraction(-this%n_, this%d_)
end function

! ----------------------------------------------------------------------
! Absolute value.
! ----------------------------------------------------------------------
impure elemental function abs_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  type(IntFraction)             :: output
  
  output = IntFraction(abs(this%n_), this%d_)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_IntFraction(this,input)
  implicit none
  
  class(IntFraction), intent(out) :: this
  type(String),       intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  
  select type(this); type is(IntFraction)
    split_string = split_line(input, '/')
    if (size(split_string)==1) then
      ! Assume the string is an integer.
      this = frac(int(split_string(1)))
    elseif (size(split_string)==2) then
      ! Assume the string is of the form 'a/b'.
      this = IntFraction(int(split_string(1)),int(split_string(2)))
    else
      call print_line('Error parsing fraction from string: '//input)
      call err()
    endif
  class default
    call err()
  end select
end subroutine

function write_IntFraction(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  type(String)                   :: output
  
  select type(this); type is(IntFraction)
    if (is_int(this)) then
     output = pad_int_to_str(this%n_)
    else
      output = pad_int_to_str(this%n_)//'/'//this%d_
    endif
  class default
    call err()
  end select
end function

impure elemental function new_IntFraction_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(IntFraction)        :: this
  
  call this%read(input)
end function
end module

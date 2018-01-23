! ======================================================================
! Integer fractions in an exact representation.
! ======================================================================
module fraction_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  type, extends(Stringable) :: IntFraction
    integer, private :: n_
    integer, private :: d_
  contains
    ! Getters
    procedure, public :: denominator
    procedure, public :: numerator
    
    ! Fraction simplification.
    procedure, private :: simplify
    
    ! Conversion from integer.
    generic,   public  :: assignment(=) => assign_IntFraction_integer
    procedure, private ::                  assign_IntFraction_integer
    
    ! Comparison.
    generic, public :: operator(==) => equality_IntFraction_IntFraction, &
                                     & equality_IntFraction_integer,     &
                                     & equality_integer_IntFraction
    procedure, private              :: equality_IntFraction_IntFraction
    procedure, private              :: equality_IntFraction_integer
    procedure, private, pass(that)  :: equality_integer_IntFraction
    
    generic, public :: operator(/=) => non_equality_IntFraction_IntFraction, &
                                     & non_equality_IntFraction_integer,     &
                                     & non_equality_integer_IntFraction
    procedure, private              :: non_equality_IntFraction_IntFraction
    procedure, private              :: non_equality_IntFraction_integer
    procedure, private, pass(that)  :: non_equality_integer_IntFraction
    
    ! Arithmetic.
    generic, public :: operator(+) => add_IntFraction_IntFraction, &
                                    & add_IntFraction_integer,     &
                                    & add_integer_IntFraction
    procedure, private             :: add_IntFraction_IntFraction
    procedure, private             :: add_IntFraction_integer
    procedure, private, pass(that) :: add_integer_IntFraction
    
    generic, public :: operator(-) => subtract_IntFraction_IntFraction, &
                                    & subtract_IntFraction_integer,     &
                                    & subtract_integer_IntFraction
    procedure, private             :: subtract_IntFraction_IntFraction
    procedure, private             :: subtract_IntFraction_integer
    procedure, private, pass(that) :: subtract_integer_IntFraction
    
    generic, public :: operator(*) => multiply_IntFraction_IntFraction, &
                                    & multiply_IntFraction_integer,     &
                                    & multiply_integer_IntFraction
    procedure, private             :: multiply_IntFraction_IntFraction
    procedure, private             :: multiply_IntFraction_integer
    procedure, private, pass(that) :: multiply_integer_IntFraction
    
    generic, public :: operator(/) => divide_IntFraction_IntFraction, &
                                    & divide_IntFraction_integer,     &
                                    & divide_integer_IntFraction
    procedure, private             :: divide_IntFraction_IntFraction
    procedure, private             :: divide_IntFraction_integer
    procedure, private, pass(that) :: divide_integer_IntFraction
    
    ! Conversions to other types.
    procedure, public  :: int => int_IntFraction
    procedure, public  :: dble => dble_IntFraction
    
    ! I/O.
    procedure, public :: str => str_IntFraction
  end type
  
  ! Constructor.
  interface IntFraction
    module procedure new_IntFraction
  end interface
  
  ! Conversions from other types.
  interface frac
    module procedure frac_character
    module procedure frac_String
    module procedure frac_integer
  end interface
  
  ! Check whether or not the fraction is an integer.
  interface is_int
    module procedure is_int_IntFraction
  end interface
  
  ! Negative.
  interface operator(-)
    module procedure negative_IntFraction
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
pure function new_IntFraction(numerator,denominator) result(this)
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
pure function numerator(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer                        :: output
  
  output = this%n_
end function

pure function denominator(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer                        :: output
  
  output = this%d_
end function

! ----------------------------------------------------------------------
! Simplifies the fraction.
! ----------------------------------------------------------------------
pure subroutine simplify(this)
  use utils_module, only : gcd
  
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
! Conversion from integer.
! ----------------------------------------------------------------------
pure subroutine assign_IntFraction_integer(output,input)
  implicit none
  
  class(IntFraction), intent(inout) :: output
  integer,            intent(in)    :: input
  
  output%n_ = input
  output%d_ = 1
end subroutine

! ----------------------------------------------------------------------
! Comparison.
! ----------------------------------------------------------------------
elemental function equality_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this%n_==that%n_ .and. this%d_==that%d_
end function

elemental function equality_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = this%n_==that .and. this%d_==1
end function

elemental function equality_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = this==that%n_ .and. that%d_==1
end function

elemental function non_equality_IntFraction_IntFraction(this,that) &
   & result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

elemental function non_equality_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

elemental function non_equality_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  logical                        :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Addition.
! ----------------------------------------------------------------------
! a/b + c/d = (ad+bc)/(bd).
elemental function add_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_+that%n_*this%d_, this%d_*that%d_)
end function

! a/b + c = (a+bc)/b.
elemental function add_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_+that*this%d_, this%d_)
end function

! a + b/c = (ac+b)/c.
elemental function add_integer_IntFraction(this,that) result(output)
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
elemental function subtract_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_-that%n_*this%d_, this%d_*that%d_)
end function

! a/b - c = (a-bc)/b.
elemental function subtract_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_-that*this%d_, this%d_)
end function

! a - b/c = (ac-b)/c.
elemental function subtract_integer_IntFraction(this,that) result(output)
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
elemental function multiply_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%n_, this%d_*that%d_)
end function

! a/b * c = ac/b.
elemental function multiply_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that, this%d_)
end function

! a * b/c = ab/c.
elemental function multiply_integer_IntFraction(this,that) result(output)
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
elemental function divide_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_, this%d_*that%n_)
end function

! a/b / c = a/(bc).
elemental function divide_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_, this%d_*that)
end function

! a / b/c = ac/b.
elemental function divide_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%d_, that%n_)
end function

! ----------------------------------------------------------------------
! Conversions from other types.
! ----------------------------------------------------------------------
elemental function frac_character(input) result(output)
  implicit none
  
  character(*), intent(in) :: input
  type(IntFraction)        :: output
  
  type(String), allocatable :: split_string(:)
  
  split_string = split(input, '/')
  output = IntFraction(int(split_string(1)),int(split_string(2)))
end function

elemental function frac_String(input) result(output)
  implicit none
  
  type(String), intent(in) :: input
  type(IntFraction)        :: output
  
  output = frac(char(input))
end function

elemental function frac_integer(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  type(IntFraction)   :: output
  
  output = input
end function

! ----------------------------------------------------------------------
! Conversions to other types.
! ----------------------------------------------------------------------
! Conversion to integer.
! As with int(real), rounds down non-integer fractions.
elemental function int_IntFraction(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer                        :: output
  
  output = this%n_/this%d_
end function

! Conversion to real(dp).
elemental function dble_IntFraction(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  real(dp)                       :: output
  
  output = real(this%n_,dp) / this%d_
end function

! ----------------------------------------------------------------------
! Whether or not the IntFraction is an integer.
! ----------------------------------------------------------------------
! Equivalent to whether or not the denominator = 1.
elemental function is_int_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  logical                       :: output
  
  output = this%d_==1
end function

! ----------------------------------------------------------------------
! Negative.
! ----------------------------------------------------------------------
elemental function negative_IntFraction(this) result(output)
  implicit none
  
  type(IntFraction), intent(in) :: this
  type(IntFraction)             :: output
  
  output = IntFraction(-this%n_, this%d_)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
pure function str_IntFraction(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  type(String)                   :: output
  
  output = this%n_//'/'//this%d_
end function
end module

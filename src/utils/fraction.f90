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
    
    ! Comparison.
    
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
    
    ! Conversion to real.
    procedure, public :: dble => dble_IntFraction
    
    ! I/O.
    procedure, public :: str => str_IntFraction
  end type
contains

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
subroutine simplify(this)
  use utils_module, only : gcd
  
  class(IntFraction), intent(inout) :: this
  
  integer :: common_denominator
  
  common_denominator = gcd(this%n_, this%d_)
  this%n_ = this%n_ / common_denominator
  this%d_ = this%d_ / common_denominator
end subroutine

! ----------------------------------------------------------------------
! Comparison.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Addition.
! ----------------------------------------------------------------------
function add_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_+that%n_*this%d_, this%d_*that%d_)
  call output%simplify()
end function

function add_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_+that*this%d_, this%d_)
end function

function add_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%d_+that%n_, that%d_)
end function

! ----------------------------------------------------------------------
! Subtraction.
! ----------------------------------------------------------------------
function subtract_IntFraction_IntFraction(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_*that%d_-that%n_*this%d_, this%d_*that%d_)
  call output%simplify()
end function

function subtract_IntFraction_integer(this,that) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  integer,            intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this%n_-that*this%d_, this%d_)
end function

function subtract_integer_IntFraction(this,that) result(output)
  implicit none
  
  integer,            intent(in) :: this
  class(IntFraction), intent(in) :: that
  type(IntFraction)              :: output
  
  output = IntFraction(this*that%d_-that%n_, that%d_)
end function

! ----------------------------------------------------------------------
! Multiplication.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Division.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Conversion to real.
! ----------------------------------------------------------------------
pure function dble_IntFraction(this) result(output)
  implicit none
  
  class(IntFraction), intent(in) :: this
  real(dp)                       :: output
  
  output = real(this%n_,dp) / this%d_
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

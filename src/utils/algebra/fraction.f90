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
  
  interface IntFraction
    ! ----------------------------------------------------------------------
    ! Constructor.
    ! ----------------------------------------------------------------------
    module function new_IntFraction(numerator,denominator) result(this) 
      integer, intent(in) :: numerator
      integer, intent(in) :: denominator
      type(IntFraction)   :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Getters.
    ! ----------------------------------------------------------------------
    impure elemental module function numerator(this) result(output) 
      class(IntFraction), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface
    impure elemental module function denominator(this) result(output) 
      class(IntFraction), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Simplifies the fraction.
    ! ----------------------------------------------------------------------
    module subroutine simplify(this) 
      class(IntFraction), intent(inout) :: this
    end subroutine
  end interface
  
  interface int
    ! ----------------------------------------------------------------------
    ! Conversions to and from other types.
    ! ----------------------------------------------------------------------
    ! Conversion to integer.
    ! As with int(real), rounds down non-integer fractions.
    impure elemental module function int_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      integer                       :: output
    end function
  end interface
  
  interface dble
    ! Conversion to real(dp).
    impure elemental module function dble_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      real(dp)                      :: output
    end function
  end interface
  
  interface frac
    ! Conversion from character(*).
    impure elemental module function frac_character(input) result(output) 
      character(*), intent(in) :: input
      type(IntFraction)        :: output
    end function
  
    ! Conversion from String.
    impure elemental module function frac_String(input) result(output) 
      type(String), intent(in) :: input
      type(IntFraction)        :: output
    end function
  
    ! Conversion from integer.
    impure elemental module function frac_integer(input) result(output) 
      integer, intent(in) :: input
      type(IntFraction)   :: output
    end function
  
    ! Conversion from numerator and denominator.
    module function frac_integers(numerator,denominator) result(output) 
      integer, intent(in) :: numerator
      integer, intent(in) :: denominator
      type(IntFraction)   :: output
    end function
  end interface
  
  interface operator(==)
    ! ----------------------------------------------------------------------
    ! Comparison.
    ! ----------------------------------------------------------------------
    impure elemental module function equality_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function equality_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function equality_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(/=)
    impure elemental module function non_equality_IntFraction_IntFraction( &
       & this,that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function non_equality_IntFraction_integer(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function non_equality_integer_IntFraction(this, &
       & that) result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(<)
    impure elemental module function lt_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function lt_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function lt_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(>)
    impure elemental module function gt_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function gt_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function gt_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(<=)
    impure elemental module function le_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function le_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function le_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(>=)
    impure elemental module function ge_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function ge_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    impure elemental module function ge_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(+)
    ! ----------------------------------------------------------------------
    ! Addition.
    ! ----------------------------------------------------------------------
    ! a/b + c/d = (ad+bc)/(bd).
    impure elemental module function add_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a/b + c = (a+bc)/b.
    impure elemental module function add_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a + b/c = (ac+b)/c.
    impure elemental module function add_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface operator(-)
    ! ----------------------------------------------------------------------
    ! Subtraction.
    ! ----------------------------------------------------------------------
    ! a/b + c/d = (ad-bc)/(bd).
    impure elemental module function subtract_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a/b - c = (a-bc)/b.
    impure elemental module function subtract_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a - b/c = (ac-b)/c.
    impure elemental module function subtract_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Multiplication.
    ! ----------------------------------------------------------------------
    ! a/b * c/d = (ac)/(bd).
    impure elemental module function multiply_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a/b * c = ac/b.
    impure elemental module function multiply_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a * b/c = ab/c.
    impure elemental module function multiply_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface operator(/)
    ! ----------------------------------------------------------------------
    ! Division.
    ! ----------------------------------------------------------------------
    ! a/b / c/d = (ad)/(bc).
    impure elemental module function divide_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a/b / c = a/(bc).
    impure elemental module function divide_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    ! a / b/c = ac/b.
    impure elemental module function divide_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface is_int
    ! ----------------------------------------------------------------------
    ! Whether or not the IntFraction is an integer.
    ! ----------------------------------------------------------------------
    ! Equivalent to whether or not the denominator = 1.
    impure elemental module function is_int_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      logical                       :: output
    end function
  end interface
  
  interface modulo
    ! ----------------------------------------------------------------------
    ! A fraction modulo an integer.
    ! ----------------------------------------------------------------------
    impure elemental module function modulo_IntFraction_integer(this,that) &
       & result(output) 
      type(IntFraction), intent(in) :: this
      integer,           intent(in) :: that
      type(IntFraction)             :: output
    end function
  end interface
  
  interface operator(-)
    ! ----------------------------------------------------------------------
    ! Negative.
    ! ----------------------------------------------------------------------
    impure elemental module function negative_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      type(IntFraction)             :: output
    end function
  end interface
  
  interface abs
    ! ----------------------------------------------------------------------
    ! Absolute value.
    ! ----------------------------------------------------------------------
    impure elemental module function abs_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      type(IntFraction)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_IntFraction(this,input) 
      class(IntFraction), intent(out) :: this
      type(String),       intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_IntFraction(this) result(output) 
      class(IntFraction), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface IntFraction
    impure elemental module function new_IntFraction_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(IntFraction)        :: this
    end function
  end interface
end module

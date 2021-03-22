!> Provides the [[IntFraction(type)]] class and related methods.
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
  
  !> Stores a fraction in exact representation.
  type, extends(Stringable) :: IntFraction
    !> The numerator.
    integer, private :: n_
    !> The denominator.
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
    ! [[IntFraction(type)]] constructor. Simplifies the fraction if possible.
    module function new_IntFraction(numerator,denominator) result(this) 
      integer, intent(in) :: numerator
      integer, intent(in) :: denominator
      type(IntFraction)   :: this
    end function
  end interface
  
  interface
    !> Returns the numerator of the (simplified) fraction.
    impure elemental module function numerator(this) result(output) 
      class(IntFraction), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface
    !> Returns the denominator of the (simplified) fraction.
    impure elemental module function denominator(this) result(output) 
      class(IntFraction), intent(in) :: this
      integer                        :: output
    end function
  end interface
  
  interface
    !> Simplifies the fraction.
    module subroutine simplify(this) 
      class(IntFraction), intent(inout) :: this
    end subroutine
  end interface
  
  interface int
    !> Conversion from [[IntFraction(type)]] to `integer`.
    !> Rounds non-integer fractions towards zero, consistent with `int(real)`.
    impure elemental module function int_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      integer                       :: output
    end function
  end interface
  
  interface dble
    !> Conversion from [[IntFraction(type)]] to `real(dp)`.
    impure elemental module function dble_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      real(dp)                      :: output
    end function
  end interface
  
  interface frac
    !> Conversion from `character(*)` to [[IntFraction(type)]].
    impure elemental module function frac_character(input) result(output) 
      character(*), intent(in) :: input
      type(IntFraction)        :: output
    end function
  
    !> Conversion from [[String(type)]] to [[IntFraction(type)]].
    impure elemental module function frac_String(input) result(output) 
      type(String), intent(in) :: input
      type(IntFraction)        :: output
    end function
  
    !> Conversion from `integer` to [[IntFraction(type)]].
    impure elemental module function frac_integer(input) result(output) 
      integer, intent(in) :: input
      type(IntFraction)   :: output
    end function
  
    !> Conversion from numerator and denominator to [[IntFraction(type)]].
    module function frac_integers(numerator,denominator) result(output) 
      integer, intent(in) :: numerator
      integer, intent(in) :: denominator
      type(IntFraction)   :: output
    end function
  end interface
  
  interface operator(==)
    !> Equality between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function equality_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    !> Equality between [[IntFraction(type)]] and `integer`.
    impure elemental module function equality_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    !> Equality between [[IntFraction(type)]] and `integer`.
    impure elemental module function equality_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(/=)
    !> Non-equality between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function non_equality_IntFraction_IntFraction( &
       & this,that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    !> Non-equality between [[IntFraction(type)]] and `integer`.
    impure elemental module function non_equality_IntFraction_integer(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    !> Non-equality between [[IntFraction(type)]] and `integer`.
    impure elemental module function non_equality_integer_IntFraction(this, &
       & that) result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(<)
    !> `<` comparison between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function lt_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    !> `<` comparison between [[IntFraction(type)]] and `integer`.
    impure elemental module function lt_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    !> `<` comparison between `integer` and [[IntFraction(type)]].
    impure elemental module function lt_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(>)
    !> `>` comparison between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function gt_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    !> `>` comparison between [[IntFraction(type)]] and `integer`.
    impure elemental module function gt_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    !> `>` comparison between `integer` and [[IntFraction(type)]].
    impure elemental module function gt_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(<=)
    !> `<=` comparison between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function le_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    !> `<=` comparison between [[IntFraction(type)]] and `integer`.
    impure elemental module function le_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    !> `<=` comparison between `integer` and [[IntFraction(type)]].
    impure elemental module function le_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(>=)
    !> `>=` comparison between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function ge_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  
    !> `>=` comparison between [[IntFraction(type)]] and `integer`.
    impure elemental module function ge_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      logical                        :: output
    end function
  
    !> `>=` comparison between `integer` and [[IntFraction(type)]].
    impure elemental module function ge_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      logical                        :: output
    end function
  end interface
  
  interface operator(+)
    !> Addition between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function add_IntFraction_IntFraction(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Addition between [[IntFraction(type)]] and `integer`.
    impure elemental module function add_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Addition between `integer` and [[IntFraction(type)]].
    impure elemental module function add_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface operator(-)
    !> Subtraction between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function subtract_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Subtraction between [[IntFraction(type)]] and `integer`.
    impure elemental module function subtract_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Subtraction between `integer` and [[IntFraction(type)]].
    impure elemental module function subtract_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface operator(*)
    !> Multiplication between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function multiply_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Multiplication between [[IntFraction(type)]] and `integer`.
    impure elemental module function multiply_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Multiplication between `integer` and [[IntFraction(type)]].
    impure elemental module function multiply_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface operator(/)
    !> Division between [[IntFraction(type)]] and [[IntFraction(type)]].
    impure elemental module function divide_IntFraction_IntFraction(this, &
       & that) result(output) 
      class(IntFraction), intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Division between [[IntFraction(type)]] and `integer`.
    impure elemental module function divide_IntFraction_integer(this,that) &
       & result(output) 
      class(IntFraction), intent(in) :: this
      integer,            intent(in) :: that
      type(IntFraction)              :: output
    end function
  
    !> Division between `integer` and [[IntFraction(type)]].
    impure elemental module function divide_integer_IntFraction(this,that) &
       & result(output) 
      integer,            intent(in) :: this
      class(IntFraction), intent(in) :: that
      type(IntFraction)              :: output
    end function
  end interface
  
  interface is_int
    !> Returns `true` if the [[IntFraction(type)]] is an integer,
    !>   and `false` otherwise.
    impure elemental module function is_int_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      logical                       :: output
    end function
  end interface
  
  interface modulo
    !> Returns `dividend` modulo `divisor`.
    !> The sign of the output matches the sign of `divisor`.
    impure elemental module function modulo_IntFraction_integer(dividend, &
       & divisor) result(output) 
      type(IntFraction), intent(in) :: dividend
      integer,           intent(in) :: divisor
      type(IntFraction)             :: output
    end function
  end interface
  
  interface operator(-)
    !> `-` negative operator for [[IntFraction(type)]].
    impure elemental module function negative_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      type(IntFraction)             :: output
    end function
  end interface
  
  interface abs
    ! Resturns the absolute value of an [[IntFraction(type)]].
    impure elemental module function abs_IntFraction(this) result(output) 
      type(IntFraction), intent(in) :: this
      type(IntFraction)             :: output
    end function
  end interface
  
  interface
    !> Converts a [[String(type)]] to [[IntFraction(type)]].
    module subroutine read_IntFraction(this,input) 
      class(IntFraction), intent(out) :: this
      type(String),       intent(in)  :: input
    end subroutine
  end interface
  
  interface
    !> Converts an [[IntFraction(type)]] to [[String(type)]].
    module function write_IntFraction(this) result(output) 
      class(IntFraction), intent(in) :: this
      type(String)                   :: output
    end function
  end interface
  
  interface IntFraction
    !> Converts a [[String(type)]] to [[IntFraction(type)]].
    impure elemental module function new_IntFraction_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(IntFraction)        :: this
    end function
  end interface
end module

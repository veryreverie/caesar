! ======================================================================
! Complex numbers with fractional real and imaginary parts.
! ======================================================================
module caesar_fraction_complex_module
  use caesar_io_module
  
  use caesar_fraction_module
  use caesar_integer_complex_module
  implicit none
  
  private
  
  public :: FractionComplex
  public :: real
  public :: aimag
  public :: conjg
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  
  type, extends(Stringable) :: FractionComplex
    type(IntFraction), private :: real_
    type(IntFraction), private :: imag_
  contains
    procedure, public :: read  => read_FractionComplex
    procedure, public :: write => write_FractionComplex
  end type
  
  interface FractionComplex
    ! Constructor.
    module function new_FractionComplex(real_part,imaginary_part) result(this) 
      type(IntFraction), intent(in)           :: real_part
      type(IntFraction), intent(in), optional :: imaginary_part
      type(FractionComplex)                   :: this
    end function
  end interface
  
  interface real
    ! Complex functions (real,aimag and conjg) .
    module function real_FractionComplex(this) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntFraction)                 :: output
    end function
  end interface
  
  interface aimag
    module function aimag_FractionComplex(this) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntFraction)                 :: output
    end function
  end interface
  
  interface conjg
    module function conjg_FractionComplex(this) result(output) 
      type(FractionComplex), intent(in) :: this
      type(FractionComplex)             :: output
    end function
  end interface
  
  interface operator(+)
    ! Arithmetic.
    impure elemental module function add_FractionComplex_FractionComplex(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function add_FractionComplex_IntComplex(this, &
       & that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntComplex),      intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function add_IntComplex_FractionComplex(this, &
       & that) result(output) 
      type(IntComplex),      intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function add_FractionComplex_IntFraction(this, &
       & that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntFraction),     intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function add_IntFraction_FractionComplex(this, &
       & that) result(output) 
      type(IntFraction),     intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function add_FractionComplex_integer(this,that) &
       & result(output) 
      type(FractionComplex), intent(in) :: this
      integer,               intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function add_integer_FractionComplex(this,that) &
       & result(output) 
      integer,               intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function subtract_FractionComplex_FractionComplex(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function subtract_FractionComplex_IntComplex(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntComplex),      intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function subtract_IntComplex_FractionComplex(this,that) result(output) 
      type(IntComplex),      intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function subtract_FractionComplex_IntFraction(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntFraction),     intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function subtract_IntFraction_FractionComplex(this,that) result(output) 
      type(IntFraction),     intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function subtract_FractionComplex_integer(this, &
       & that) result(output) 
      type(FractionComplex), intent(in) :: this
      integer,               intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function subtract_integer_FractionComplex(this, &
       & that) result(output) 
      integer,               intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  end interface
  
  interface operator(*)
    impure elemental module function multiply_FractionComplex_FractionComplex(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function multiply_FractionComplex_IntComplex(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntComplex),      intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function multiply_IntComplex_FractionComplex(this,that) result(output) 
      type(IntComplex),      intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function multiply_FractionComplex_IntFraction(this,that) result(output) 
      type(FractionComplex), intent(in) :: this
      type(IntFraction),     intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function multiply_IntFraction_FractionComplex(this,that) result(output) 
      type(IntFraction),     intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function multiply_FractionComplex_integer(this, &
       & that) result(output) 
      type(FractionComplex), intent(in) :: this
      integer,               intent(in) :: that
      type(FractionComplex)             :: output
    end function
  
    impure elemental module function multiply_integer_FractionComplex(this, &
       & that) result(output) 
      integer,               intent(in) :: this
      type(FractionComplex), intent(in) :: that
      type(FractionComplex)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_FractionComplex(this,input) 
      class(FractionComplex), intent(out) :: this
      type(String),           intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_FractionComplex(this) result(output) 
      class(FractionComplex), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface FractionComplex
    impure elemental module function new_FractionComplex_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(FractionComplex)    :: this
    end function
  end interface
end module

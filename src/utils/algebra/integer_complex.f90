! ======================================================================
! Complex numbers with integer real and imaginary parts.
! ======================================================================
module caesar_integer_complex_module
  use caesar_io_module
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
    procedure, public :: read  => read_IntComplex
    procedure, public :: write => write_IntComplex
  end type
  
  interface IntComplex
    ! Constructor.
    module function new_IntComplex(real_part,imaginary_part) result(this) 
      integer, intent(in)           :: real_part
      integer, intent(in), optional :: imaginary_part
      type(IntComplex)              :: this
    end function
  end interface
  
  interface real
    ! Complex functions (real,aimag and conjg) .
    module function real_IntComplex(this) result(output) 
      type(IntComplex), intent(in) :: this
      integer                      :: output
    end function
  end interface
  
  interface aimag
    module function aimag_IntComplex(this) result(output) 
      type(IntComplex), intent(in) :: this
      integer                      :: output
    end function
  end interface
  
  interface conjg
    module function conjg_IntComplex(this) result(output) 
      type(IntComplex), intent(in) :: this
      type(IntComplex)             :: output
    end function
  end interface
  
  interface operator(+)
    ! Arithmetic.
    impure elemental module function add_IntComplex_IntComplex(this,that) &
       & result(output) 
      type(IntComplex), intent(in) :: this
      type(IntComplex), intent(in) :: that
      type(IntComplex)             :: output
    end function
  
    impure elemental module function add_IntComplex_integer(this,that) &
       & result(output) 
      type(IntComplex), intent(in) :: this
      integer,          intent(in) :: that
      type(IntComplex)             :: output
    end function
  
    impure elemental module function add_integer_IntComplex(this,that) &
       & result(output) 
      integer,          intent(in) :: this
      type(IntComplex), intent(in) :: that
      type(IntComplex)             :: output
    end function
  end interface
  
  interface operator(-)
    impure elemental module function negative_IntComplex(this) result(output)
      type(IntComplex), intent(in) :: this
      type(IntComplex)             :: output
    end function
    
    impure elemental module function subtract_IntComplex_IntComplex(this, &
       & that) result(output) 
      type(IntComplex), intent(in) :: this
      type(IntComplex), intent(in) :: that
      type(IntComplex)             :: output
    end function
  
    impure elemental module function subtract_IntComplex_integer(this,that) &
       & result(output) 
      type(IntComplex), intent(in) :: this
      integer,          intent(in) :: that
      type(IntComplex)             :: output
    end function
  
    impure elemental module function subtract_integer_IntComplex(this,that) &
       & result(output) 
      integer,          intent(in) :: this
      type(IntComplex), intent(in) :: that
      type(IntComplex)             :: output
    end function
  end interface
  
  interface operator(*)
    impure elemental module function multiply_IntComplex_IntComplex(this, &
       & that) result(output) 
      type(IntComplex), intent(in) :: this
      type(IntComplex), intent(in) :: that
      type(IntComplex)             :: output
    end function
  
    impure elemental module function multiply_IntComplex_integer(this,that) &
       & result(output) 
      type(IntComplex), intent(in) :: this
      integer,          intent(in) :: that
      type(IntComplex)             :: output
    end function
  
    impure elemental module function multiply_integer_IntComplex(this,that) &
       & result(output) 
      integer,          intent(in) :: this
      type(IntComplex), intent(in) :: that
      type(IntComplex)             :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_IntComplex(this,input) 
      class(IntComplex), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_IntComplex(this) result(output) 
      class(IntComplex), intent(in) :: this
      type(String)                  :: output
    end function
  end interface
  
  interface IntComplex
    impure elemental module function new_IntComplex_String(input) result(this) 
      type(String), intent(in) :: input
      type(IntComplex)         :: this
    end function
  end interface
end module

! ======================================================================
! Complex numbers with fractional real and imaginary parts.
! ======================================================================
module fraction_complex_submodule
  use io_module
  
  use fraction_submodule
  use integer_complex_submodule
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
    module procedure new_FractionComplex
    module procedure new_FractionComplex_String
  end interface
  
  interface real
    module procedure real_FractionComplex
  end interface
  
  interface aimag
    module procedure aimag_FractionComplex
  end interface
  
  interface conjg
    module procedure conjg_FractionComplex
  end interface
  
  interface operator(+)
    module procedure add_FractionComplex_FractionComplex
    module procedure add_FractionComplex_IntComplex
    module procedure add_IntComplex_FractionComplex
    module procedure add_FractionComplex_IntFraction
    module procedure add_IntFraction_FractionComplex
    module procedure add_FractionComplex_integer
    module procedure add_integer_FractionComplex
  end interface
  
  interface operator(-)
    module procedure subtract_FractionComplex_FractionComplex
    module procedure subtract_FractionComplex_IntComplex
    module procedure subtract_IntComplex_FractionComplex
    module procedure subtract_FractionComplex_IntFraction
    module procedure subtract_IntFraction_FractionComplex
    module procedure subtract_FractionComplex_integer
    module procedure subtract_integer_FractionComplex
  end interface
  
  interface operator(*)
    module procedure multiply_FractionComplex_FractionComplex
    module procedure multiply_FractionComplex_IntComplex
    module procedure multiply_IntComplex_FractionComplex
    module procedure multiply_FractionComplex_IntFraction
    module procedure multiply_IntFraction_FractionComplex
    module procedure multiply_FractionComplex_integer
    module procedure multiply_integer_FractionComplex
  end interface
contains

! Constructor.
function new_FractionComplex(real_part,imaginary_part) result(this)
  implicit none
  
  type(IntFraction), intent(in)           :: real_part
  type(IntFraction), intent(in), optional :: imaginary_part
  type(FractionComplex)                   :: this
  
  this%real_ = real_part
  if (present(imaginary_part)) then
    this%imag_ = imaginary_part
  else
    this%imag_ = 0
  endif
end function

! Complex functions (real, aimag and conjg).
function real_FractionComplex(this) result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntFraction)                 :: output
  
  output = this%real_
end function

function aimag_FractionComplex(this) result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntFraction)                 :: output
  
  output = this%imag_
end function

function conjg_FractionComplex(this) result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(FractionComplex)             :: output
  
  output = FractionComplex(real(this),-aimag(this))
end function


! Arithmetic.
impure elemental function add_FractionComplex_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end function

impure elemental function add_FractionComplex_IntComplex(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntComplex),      intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end function

impure elemental function add_IntComplex_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntComplex),      intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) + real(that)
  output%imag_ = aimag(this) + aimag(that)
end function

impure elemental function add_FractionComplex_IntFraction(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntFraction),     intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) + that
  output%imag_ = aimag(this)
end function

impure elemental function add_IntFraction_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntFraction),     intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = this + real(that)
  output%imag_ =        aimag(that)
end function

impure elemental function add_FractionComplex_integer(this,that) result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  integer,               intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) + that
  output%imag_ = aimag(this)
end function

impure elemental function add_integer_FractionComplex(this,that) result(output)
  implicit none
  
  integer,               intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = this + real(that)
  output%imag_ =        aimag(that)
end function

impure elemental function subtract_FractionComplex_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end function

impure elemental function subtract_FractionComplex_IntComplex(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntComplex),      intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end function

impure elemental function subtract_IntComplex_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntComplex),      intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) - real(that)
  output%imag_ = aimag(this) - aimag(that)
end function

impure elemental function subtract_FractionComplex_IntFraction(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntFraction),     intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) - that
  output%imag_ = aimag(this)
end function

impure elemental function subtract_IntFraction_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntFraction),     intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = this - real(that)
  output%imag_ =      - aimag(that)
end function

impure elemental function subtract_FractionComplex_integer(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  integer,               intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) - that
  output%imag_ = aimag(this)
end function

impure elemental function subtract_integer_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  integer,               intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = this - real(that)
  output%imag_ =      - aimag(that)
end function

impure elemental function multiply_FractionComplex_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end function

impure elemental function multiply_FractionComplex_IntComplex(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntComplex),      intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end function

impure elemental function multiply_IntComplex_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntComplex),      intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this)*real(that) - aimag(this)*aimag(that)
  output%imag_ = real(this)*aimag(that) + aimag(this)*real(that)
end function

impure elemental function multiply_FractionComplex_IntFraction(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  type(IntFraction),     intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) * that
  output%imag_ = aimag(this) * that
end function

impure elemental function multiply_IntFraction_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  type(IntFraction),     intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = this * real(that)
  output%imag_ = this * aimag(that)
end function

impure elemental function multiply_FractionComplex_integer(this,that) &
   & result(output)
  implicit none
  
  type(FractionComplex), intent(in) :: this
  integer,               intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = real(this) * that
  output%imag_ = aimag(this) * that
end function

impure elemental function multiply_integer_FractionComplex(this,that) &
   & result(output)
  implicit none
  
  integer,               intent(in) :: this
  type(FractionComplex), intent(in) :: that
  type(FractionComplex)             :: output
  
  output%real_ = this * real(that)
  output%imag_ = this * aimag(that)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_FractionComplex(this,input)
  implicit none
  
  class(FractionComplex), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  type(IntFraction)         :: real_part
  type(IntFraction)         :: imag_part
  
  select type(this); type is(FractionComplex)
    split_string = split_line(input,'+')
    if (size(split_string)==2) then
      ! input is of the form "a+bi" or "-a+bi".
      real_part = frac(split_string(1))
      imag_part = frac(slice(split_string(2),1,len(split_string(2))-1))
    elseif (size(split_string)==1) then
      ! input is of the form "a-bi" or "-a-bi".
      split_string = split_line(input,'-')
      if (size(split_string)==2) then
        if (slice(input,1,1)=='-') then
          real_part = -frac(split_string(1))
        else
          real_part = frac(split_string(1))
        endif
        imag_part = frac(slice(split_string(2),1,len(split_string(2))-1))
      else
        call print_line(ERROR//': Unable to read IntComplex from string: '// &
           & input)
        call err()
      endif
    else
      call print_line(ERROR//': Unable to read IntComplex from string: '// &
         & input)
      call err()
    endif
    
    this = FractionComplex(real_part,imag_part)
  end select
end subroutine

function write_FractionComplex(this) result(output)
  implicit none
  
  class(FractionComplex), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(FractionComplex)
    ! N.B. abs() is called in both cases to stop +/-0 from breaking formatting.
    if (aimag(this)>=0) then
      output = real(this)//'+'//abs(aimag(this))//'i'
    else
      output = real(this)//'-'//abs(aimag(this))//'i'
    endif
  end select
end function

impure elemental function new_FractionComplex_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(FractionComplex)    :: this
  
  this = input
end function
end module

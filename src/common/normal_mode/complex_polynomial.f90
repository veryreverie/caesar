! ======================================================================
! The building blocks of basis functions in complex co-ordinates.
! ======================================================================
module complex_polynomial_submodule
  use utils_module
  
  use complex_mode_submodule
  use complex_single_mode_displacement_submodule
  use complex_mode_displacement_submodule
  implicit none
  
  private
  
  public :: ComplexUnivariate
  public :: ComplexMonomial
  public :: ComplexPolynomial
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: conjg
  
  ! --------------------------------------------------
  ! Types, and conversions between types.
  ! --------------------------------------------------
  
  ! It is desirable to be able to convert:
  !    univariates -> monomials -> polynomials
  ! This is acheived by extending from the classes ComplexPolynomialable and
  !    ComplexMonomialable, representing types which can be converted to
  !    ComplexPolynomial and ComplexMonomial respectively.
  
  type, abstract, extends(Stringable) :: ComplexPolynomialable
  contains
    procedure(to_ComplexPolynomial_ComplexPolynomialable), deferred, public &
       & :: to_ComplexPolynomial
  end type
  
  type, abstract, extends(ComplexPolynomialable) :: ComplexMonomialable
  contains
    procedure(to_ComplexMonomial_ComplexMonomialable), deferred, public :: &
       & to_ComplexMonomial
  end type
  
  type, extends(ComplexMonomialable) :: ComplexUnivariate
    integer :: id
    integer :: paired_id
    integer :: power
  contains
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexUnivariate
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexUnivariate
    
    procedure, public :: evaluate   => evaluate_ComplexUnivariate
    
    procedure, public :: read  => read_ComplexUnivariate
    procedure, public :: write => write_ComplexUnivariate
  end type
  
  interface ComplexUnivariate
    module procedure new_ComplexUnivariate
  end interface
  
  type, extends(ComplexMonomialable) :: ComplexMonomial
    complex(dp)                          :: coefficient
    type(ComplexUnivariate), allocatable :: modes(:)
  contains
    procedure, public :: to_ComplexMonomial   => &
       & to_ComplexMonomial_ComplexMonomial
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexMonomial
    
    procedure, public :: evaluate   => evaluate_ComplexMonomial
    procedure, public :: derivative => derivative_ComplexMonomial
    
    procedure, public :: read  => read_ComplexMonomial
    procedure, public :: write => write_ComplexMonomial
  end type
  
  interface ComplexMonomial
    module procedure new_ComplexMonomial
  end interface
  
  type, extends(ComplexPolynomialable) :: ComplexPolynomial
    type(ComplexMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_ComplexPolynomial => &
       & to_ComplexPolynomial_ComplexPolynomial
    
    procedure, public :: evaluate   => evaluate_ComplexPolynomial
    procedure, public :: derivative => derivative_ComplexPolynomial
    
    procedure, public :: read  => read_ComplexPolynomial
    procedure, public :: write => write_ComplexPolynomial
  end type
  
  interface ComplexPolynomial
    module procedure new_ComplexPolynomial
  end interface
  
  abstract interface
    function to_ComplexPolynomial_ComplexPolynomialable(this) result(output)
      import ComplexPolynomial
      import ComplexPolynomialable
      implicit none
      
      class(ComplexPolynomialable), intent(in) :: this
      type(ComplexPolynomial)                  :: output
    end function
    
    function to_ComplexMonomial_ComplexMonomialable(this) result(output)
      import ComplexMonomial
      import ComplexMonomialable
      implicit none
      
      class(ComplexMonomialable), intent(in) :: this
      type(ComplexMonomial)                  :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Operations involving types.
  ! --------------------------------------------------
  
  interface size
    module procedure size_ComplexMonomial
    module procedure size_ComplexPolynomial
  end interface
  
  interface conjg
    module procedure conjg_ComplexUnivariate
    module procedure conjg_ComplexMonomial
    module procedure conjg_ComplexPolynomial
  end interface
  
  interface operator(*)
    module procedure multiply_ComplexMonomial_real
    module procedure multiply_real_ComplexMonomial
    module procedure multiply_ComplexMonomial_complex
    module procedure multiply_complex_ComplexMonomial
    
    module procedure multiply_ComplexMonomialable_ComplexMonomialable
  end interface
  
  interface operator(/)
    module procedure divide_ComplexMonomial_real
    module procedure divide_ComplexMonomial_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexPolynomialable_ComplexPolynomialable
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_ComplexUnivariate(id,paired_id,power) result(this)
  implicit none
  
  integer, intent(in)     :: id
  integer, intent(in)     :: paired_id
  integer, intent(in)     :: power
  type(ComplexUnivariate) :: this
  
  this%id        = id
  this%paired_id = paired_id
  this%power     = power
end function

function new_ComplexMonomial(coefficient,modes) result(this)
  implicit none
  
  complex(dp),             intent(in) :: coefficient
  type(ComplexUnivariate), intent(in) :: modes(:)
  type(ComplexMonomial)               :: this
  
  this%coefficient = coefficient
  this%modes       = modes
end function

function new_ComplexPolynomial(terms) result(this)
  implicit none
  
  type(ComplexMonomial), intent(in) :: terms(:)
  type(ComplexPolynomial)           :: this
  
  this%terms = terms
end function

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
function to_ComplexMonomial_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexMonomial)                :: output
  
  output = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                          & modes       = [this])
end function

function to_ComplexPolynomial_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  output = ComplexPolynomial([this%to_ComplexMonomial()])
end function

function to_ComplexMonomial_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)              :: output
  
  output = this
end function

function to_ComplexPolynomial_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(ComplexPolynomial)            :: output
  
  output = ComplexPolynomial([this])
end function

function to_ComplexPolynomial_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)              :: output
  
  output = this
end function

! ----------------------------------------------------------------------
! Operations involving types.
! ----------------------------------------------------------------------

! The number of modes in a monomial, or terms in a polynomial.
function size_ComplexMonomial(this) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  integer                           :: output
  
  output = size(this%modes)
end function

function size_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  integer                              :: output
  
  output = size(this%terms)
end function

! Find the conjugate of a univariate or monomial.
impure elemental function conjg_ComplexUnivariate(this) result(output)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: this
  type(ComplexUnivariate)             :: output
  
  output = ComplexUnivariate( id        = this%paired_id, &
                            & paired_id = this%id,        &
                            & power     = this%power)
end function

impure elemental function conjg_ComplexMonomial(this) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  type(ComplexMonomial)             :: output
  
  output = ComplexMonomial( coefficient = conjg(this%coefficient), &
                          & modes       = conjg(this%modes))
  output%modes = output%modes(sort(output%modes%id))
end function

impure elemental function conjg_ComplexPolynomial(this) result(output)
  implicit none
  
  type(ComplexPolynomial), intent(in) :: this
  type(ComplexPolynomial)             :: output
  
  output = ComplexPolynomial(conjg(this%terms))
end function

! Evaluate a univariate, monomial or polynomial at a given displacement.
function evaluate_ComplexUnivariate(this,displacement) result(output)
  implicit none
  
  class(ComplexUnivariate),            intent(in) :: this
  type(ComplexSingleModeDisplacement), intent(in) :: displacement
  complex(dp)                                     :: output
  
  if (this%id/=displacement%id) then
    call print_line(CODE_ERROR//': Trying to evaluate a univariate at an &
       &incompatible displacement.')
    call err()
  endif
  
  output = displacement%displacement**this%power
end function

function evaluate_ComplexMonomial(this,displacement) result(output)
  implicit none
  
  class(ComplexMonomial),        intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  integer :: i,j
  
  output = this%coefficient
  
  do i=1,size(this)
    ! Find the mode in the displacement which matches that in the monomial.
    j = first(displacement%displacements%id==this%modes(i)%id,default=0)
    
    ! If the mode is not present in the displacement, then the displacement
    !    is zero. As such, the monomial is zero. (0**n=0 if n>0).
    if (j==0) then
      output = 0.0_dp
      return
    endif
    
    ! If the mode is present in both, evaluate the univariate at the
    !    displacement.
    output = output * this%modes(i)%evaluate(displacement%displacements(j))
  enddo
end function

function evaluate_ComplexPolynomial(this,displacement) result(output)
  implicit none
  
  class(ComplexPolynomial),      intent(in) :: this
  type(ComplexModeDisplacement), intent(in) :: displacement
  complex(dp)                               :: output
  
  integer :: i
  
  output = 0
  
  do i=1,size(this)
    output = output + this%terms(i)%evaluate(displacement)
  enddo
end function

! Take the derivative of a monomial or polynomial in the direction of the
!    given mode.
function derivative_ComplexMonomial(this,mode_id) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  integer,                intent(in) :: mode_id
  type(ComplexMonomial)              :: output
  
  integer :: i
  
  ! Find the univariate corresponding to the mode.
  i = first(this%modes%id==mode_id,default=0)
  
  if (i==0) then
    ! If the mode is not present, then the derivative is zero.
    output%coefficient = 0
    output%modes = [ComplexUnivariate::]
  else
    ! If the mode is present, then u^n -> n*u^(n-1).
    output = this
    output%coefficient = output%coefficient * output%modes(i)%power
    output%modes(i)%power = output%modes(i)%power - 1
    
    ! If n-1=0, remove that univariate.
    if (output%modes(i)%power==0) then
      output%modes = [output%modes(:i-1), output%modes(i+1:)]
    endif
  endif
end function

function derivative_ComplexPolynomial(this,mode_id) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  integer,                  intent(in) :: mode_id
  type(ComplexPolynomial)              :: output
  
  integer :: i,ialloc
  
  allocate(output%terms(size(this)), stat=ialloc); call err(ialloc)
  
  ! Take derivatives, term by term.
  do i=1,size(this)
    output%terms(i) = this%terms(i)%derivative(mode_id)
  enddo
  
  ! Remove the terms which are now zero.
  output%terms = output%terms(filter(output%terms,non_zero))
contains
  ! A Lambda for identifying non-zero monomials.
  function non_zero(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(ComplexMonomial)
      output = size(input)/=0
    end select
  end function
end function

! Multiplication and division by scalars.
impure elemental function multiply_ComplexMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  real(dp),              intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient * that
end function

impure elemental function multiply_real_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = that
  output%coefficient = output%coefficient * this
end function

impure elemental function multiply_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient * that
end function

impure elemental function multiply_complex_ComplexMonomial(this,that) &
   & result(output)
  implicit none
  
  complex(dp),           intent(in) :: this
  type(ComplexMonomial), intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = that
  output%coefficient = output%coefficient * this
end function

impure elemental function divide_ComplexMonomial_real(this,that) result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  real(dp),              intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient / that
end function

impure elemental function divide_ComplexMonomial_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexMonomial), intent(in) :: this
  complex(dp),           intent(in) :: that
  type(ComplexMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient / that
end function

! Multiplication between Monomials and Monomial-like types.
function multiply_ComplexMonomialable_ComplexMonomialable(this,that) &
   & result(output)
  implicit none
  
  class(ComplexMonomialable), intent(in) :: this
  class(ComplexMonomialable), intent(in) :: that
  type(ComplexMonomial)                  :: output
  
  type(ComplexMonomial) :: this_monomial
  type(ComplexMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_ComplexMonomial()
  that_monomial = that%to_ComplexMonomial()
  
  output%coefficient = this_monomial%coefficient * that_monomial%coefficient
  
  if (size(this_monomial)==0) then
    output%modes = that_monomial%modes
  elseif (size(that_monomial)==0) then
    output%modes = this_monomial%modes
  else
    i_this = 1
    i_that = 1
    i_out = 0
    allocate( output%modes(size(this_monomial)+size(that_monomial)), &
            & stat=ialloc); call err(ialloc)
    do while(i_this<=size(this_monomial) .and. i_that<=size(that_monomial))
      i_out = i_out + 1
      if (i_this>size(this_monomial)) then
        output%modes(i_out) = that_monomial%modes(i_that)
        i_that = i_that + 1
      elseif (i_that>size(that_monomial)) then
        output%modes(i_out) = this_monomial%modes(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes(i_this)%id == &
             & that_monomial%modes(i_that)%id) then
        output%modes(i_out)%id = this_monomial%modes(i_this)%id
        output%modes(i_out)%power = this_monomial%modes(i_this)%power &
                                & + that_monomial%modes(i_that)%power
        i_this = i_this + 1
        i_that = i_that + 1
      elseif ( this_monomial%modes(i_this)%id < &
             & that_monomial%modes(i_that)%id) then
        output%modes(i_out) = this_monomial%modes(i_this)
        i_this = i_this + 1
      elseif ( this_monomial%modes(i_this)%id > &
             & that_monomial%modes(i_that)%id) then
        output%modes(i_out) = that_monomial%modes(i_that)
        i_that = i_that + 1
      else
        call err()
      endif
    enddo
    output%modes = output%modes(:i_out)
  endif
end function

! Addition between polynomials and polynomial-like types.
function add_ComplexPolynomialable_ComplexPolynomialable(this,that) &
   & result(output)
  implicit none
  
  class(ComplexPolynomialable), intent(in) :: this
  class(ComplexPolynomialable), intent(in) :: that
  type(ComplexPolynomial)                  :: output
  
  type(ComplexPolynomial) :: this_polynomial
  type(ComplexPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_ComplexPolynomial()
  that_polynomial = that%to_ComplexPolynomial()
  
  allocate( output%terms(size(this_polynomial)+size(that_polynomial)), &
          & stat=ialloc); call err(ialloc)
  output%terms(:size(this_polynomial)) = this_polynomial%terms
  no_terms = size(this_polynomial)
  do i=1,size(that_polynomial)
    j = first( this_polynomial%terms,    &
             & compare_monomial_modes,   &
             & that_polynomial%terms(j), &
             & default=0)
    if (j==0) then
      no_terms = no_terms + 1
      output%terms(no_terms) = that_polynomial%terms(i)
    else
      output%terms(j)%coefficient = output%terms(j)%coefficient &
                                & + that_polynomial%terms(i)%coefficient
    endif
  enddo
  output%terms = output%terms(:no_terms)
contains
  ! Lambda for comparing monomials.
  function compare_monomial_modes(this,that) result(output)
    implicit none
    
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(ComplexMonomial)
      select type(that); type is(ComplexMonomial)
        if (size(this%modes)/=size(that%modes)) then
          output = .false.
        else
          output = all( this%modes%id==that%modes%id .and. &
                      & this%modes%power==that%modes%power)
        endif
      end select
    end select
  end function
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexUnivariate(this,input)
  implicit none
  
  class(ComplexUnivariate), intent(out) :: this
  type(String),             intent(in) :: input
  
  type(String), allocatable :: split_string(:)
  
  select type(this); type is(ComplexUnivariate)
    split_string = split_line(input,delimiter='^')
    if (size(split_string)/=2) then
      call print_line(ERROR//': Unable to convert string to univariate.')
      call err()
    endif
    
    this%id = int(slice(split_string(1),2,len(split_string(1))))
    this%power = int(split_string(2))
  end select
end subroutine

function write_ComplexUnivariate(this) result(output)
  implicit none
  
  class(ComplexUnivariate), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ComplexUnivariate)
    output = 'u'//this%id//'^'//this%power
  end select
end function

subroutine read_ComplexMonomial(this,input)
  implicit none
  
  class(ComplexMonomial), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexMonomial)
    split_string = split_line(input,delimiter='*')
    this%coefficient = dble(split_string(1))
    allocate(this%modes(size(split_string)-1), stat=ialloc); call err(ialloc)
    do i=1,size(this%modes)
      this%modes(i) = split_string(i+1)
    enddo
  end select
end subroutine

function write_ComplexMonomial(this) result(output)
  implicit none
  
  class(ComplexMonomial), intent(in) :: this
  type(String)                       :: output
  
  integer :: i
  
  select type(this); type is(ComplexMonomial)
    output = this%coefficient
    do i=1,size(this%modes)
      output = output//'*'//this%modes(i)
    enddo
  end select
end function

subroutine read_ComplexPolynomial(this,input)
  implicit none
  
  class(ComplexPolynomial), intent(out) :: this
  type(String),             intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexPolynomial)
    line = split_line(input)
    line = line(filter(line/='+'))
    allocate(this%terms(size(line)), stat=ialloc); call err(ialloc)
    do i=1,size(line)
      this%terms(i) = line(i)
    enddo
  end select
end subroutine

function write_ComplexPolynomial(this) result(output)
  implicit none
  
  class(ComplexPolynomial), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(ComplexPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end function
end module

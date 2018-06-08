! ======================================================================
! The building blocks of basis functions in real co-ordinates.
! ======================================================================
module real_polynomial_submodule
  use utils_module
  
  use real_single_mode_vector_submodule
  use real_mode_displacement_submodule
  implicit none
  
  private
  
  public :: RealUnivariate
  public :: RealMonomial
  public :: RealPolynomial
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  
  ! --------------------------------------------------
  ! Types, and conversions between types.
  ! --------------------------------------------------
  
  ! It is desirable to be able to convert:
  !    univariates -> monomials -> polynomials
  ! This is acheived by extending from the classes ComplexPolynomialable and
  !    ComplexMonomialable, representing types which can be converted to
  !    ComplexPolynomial and ComplexMonomial respectively.
  
  type, abstract, extends(Stringable) :: RealPolynomialable
  contains
    procedure(to_RealPolynomial_RealPolynomialable), deferred, public :: &
       & to_RealPolynomial
  end type
  
  type, abstract, extends(RealPolynomialable) :: RealMonomialable
  contains
    procedure(to_RealMonomial_RealMonomialable), deferred, public :: &
       & to_RealMonomial
  end type
  
  type, extends(RealMonomialable) :: RealUnivariate
    integer          :: id
    integer, private :: paired_id_
    integer          :: power
  contains
    procedure, public :: paired_id
    
    procedure, public :: to_RealMonomial   => to_RealMonomial_RealUnivariate
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealUnivariate
    
    procedure, public :: evaluate => evaluate_RealUnivariate
    
    procedure, public :: read  => read_RealUnivariate
    procedure, public :: write => write_RealUnivariate
  end type
  
  interface RealUnivariate
    module procedure new_RealUnivariate
    module procedure new_RealUnivariate_String
  end interface
  
  type, extends(RealMonomialable) :: RealMonomial
    real(dp)                          :: coefficient
    type(RealUnivariate), allocatable :: modes(:)
  contains
    procedure, public :: to_RealMonomial   => to_RealMonomial_RealMonomial
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealMonomial
    
    procedure, public :: evaluate   => evaluate_RealMonomial
    procedure, public :: derivative => derivative_RealMonomial
    
    procedure, public :: read  => read_RealMonomial
    procedure, public :: write => write_RealMonomial
  end type
  
  interface RealMonomial
    module procedure new_RealMonomial
    module procedure new_RealMonomial_String
  end interface
  
  type, extends(RealPolynomialable) :: RealPolynomial
    type(RealMonomial), allocatable :: terms(:)
  contains
    procedure, public :: to_RealPolynomial => to_RealPolynomial_RealPolynomial
    
    procedure, public :: evaluate   => evaluate_RealPolynomial
    procedure, public :: derivative => derivative_RealPolynomial
    
    procedure, public :: read  => read_RealPolynomial
    procedure, public :: write => write_RealPolynomial
  end type
  
  interface RealPolynomial
    module procedure new_RealPolynomial
    module procedure new_RealPolynomial_String
  end interface
  
  abstract interface
    function to_RealPolynomial_RealPolynomialable(this) result(output)
      import RealPolynomial
      import RealPolynomialable
      implicit none
      
      class(RealPolynomialable), intent(in) :: this
      type(RealPolynomial)                  :: output
    end function
    
    function to_RealMonomial_RealMonomialable(this) result(output)
      import RealMonomial
      import RealMonomialable
      implicit none
      
      class(RealMonomialable), intent(in) :: this
      type(RealMonomial)                  :: output
    end function
  end interface
  
  ! --------------------------------------------------
  ! Operations involving types.
  ! --------------------------------------------------
  
  interface size
    module procedure size_RealMonomial
    module procedure size_RealPolynomial
  end interface
  
  interface operator(*)
    module procedure multiply_RealMonomial_real
    module procedure multiply_real_RealMonomial
    
    module procedure multiply_RealMonomialable_RealMonomialable
  end interface
  
  interface operator(/)
    module procedure divide_RealMonomial_real
  end interface
  
  interface operator(+)
    module procedure add_RealPolynomialable_RealPolynomialable
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_RealUnivariate(id,paired_id,power) result(this)
  implicit none
  
  integer, intent(in)           :: id
  integer, intent(in), optional :: paired_id
  integer, intent(in)           :: power
  type(RealUnivariate)          :: this
  
  this%id = id
  if (present(paired_id)) then
    this%paired_id_ = paired_id
  else
    this%paired_id_ = 0
  endif
  this%power = power
end function

function new_RealMonomial(coefficient,modes) result(this)
  implicit none
  
  real(dp),             intent(in) :: coefficient
  type(RealUnivariate), intent(in) :: modes(:)
  type(RealMonomial)               :: this
  
  this%coefficient = coefficient
  this%modes       = modes
end function

function new_RealPolynomial(terms) result(this)
  implicit none
  
  type(RealMonomial), intent(in) :: terms(:)
  type(RealPolynomial)           :: this
  
  this%terms = terms
end function

! ----------------------------------------------------------------------
! Getter for paired_id, with error checking.
! ----------------------------------------------------------------------
impure elemental function paired_id(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  integer                           :: output
  
  if (this%paired_id_==0) then
    call print_line(CODE_ERROR//': Trying to use paired_id before it has &
       & been set.')
    call err()
  endif
  
  output = this%paired_id_
end function

! ----------------------------------------------------------------------
! Conversions between types.
! ----------------------------------------------------------------------
function to_RealMonomial_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  type(RealMonomial)                :: output
  
  output = RealMonomial( coefficient = 1.0_dp, &
                       & modes       = [this])
end function

function to_RealPolynomial_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  type(RealPolynomial)              :: output
  
  output = RealPolynomial([this%to_RealMonomial()])
end function

function to_RealMonomial_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  type(RealMonomial)              :: output
  
  output = this
end function

function to_RealPolynomial_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  type(RealPolynomial)            :: output
  
  output = RealPolynomial([this])
end function

function to_RealPolynomial_RealPolynomial(this) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  type(RealPolynomial)              :: output
  
  output = this
end function

! ----------------------------------------------------------------------
! Operations involving types.
! ----------------------------------------------------------------------

! The number of modes in a monomial, or terms in a polynomial.
function size_RealMonomial(this) result(output)
  implicit none
  
  type(RealMonomial), intent(in) :: this
  integer                        :: output
  
  output = size(this%modes)
end function

function size_RealPolynomial(this) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  integer                           :: output
  
  output = size(this%terms)
end function

! Evaluate a univariate, monomial or polynomial at a given displacement.
function evaluate_RealUnivariate(this,displacement) result(output)
  implicit none
  
  class(RealUnivariate),      intent(in) :: this
  type(RealSingleModeVector), intent(in) :: displacement
  real(dp)                               :: output
  
  if (this%id/=displacement%id) then
    call print_line(CODE_ERROR//': Trying to evaluate a univariate at an &
       &incompatible displacement.')
    call err()
  endif
  
  output = displacement%magnitude**this%power
end function

function evaluate_RealMonomial(this,displacement) result(output)
  implicit none
  
  class(RealMonomial),        intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  integer :: i,j
  
  output = this%coefficient
  
  do i=1,size(this)
    ! Find the mode in the displacement which matches that in the monomial.
    j = first(displacement%vectors%id==this%modes(i)%id,default=0)
    
    ! If the mode is not present in the displacement, then the displacement
    !    is zero. As such, the monomial is zero. (0**n=0 if n>0).
    if (j==0) then
      output = 0.0_dp
      return
    endif
    
    ! If the mode is present in both, evaluate the univariate at the
    !    displacement.
    output = output * this%modes(i)%evaluate(displacement%vectors(j))
  enddo
end function

function evaluate_RealPolynomial(this,displacement) result(output)
  implicit none
  
  class(RealPolynomial),      intent(in) :: this
  type(RealModeDisplacement), intent(in) :: displacement
  real(dp)                               :: output
  
  integer :: i
  
  output = 0
  
  do i=1,size(this)
    output = output + this%terms(i)%evaluate(displacement)
  enddo
end function

! Take the derivative of a monomial or polynomial in the direction of the
!    given mode.
function derivative_RealMonomial(this,mode_id) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  integer,             intent(in) :: mode_id
  type(RealMonomial)              :: output
  
  integer :: i
  
  ! Find the univariate corresponding to the mode.
  i = first(this%modes%id==mode_id,default=0)
  
  if (i==0) then
    ! If the mode is not present, then the derivative is zero.
    output%coefficient = 0
    output%modes = [RealUnivariate::]
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

function derivative_RealPolynomial(this,mode_id) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  integer,               intent(in) :: mode_id
  type(RealPolynomial)              :: output
  
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
    
    select type(input); type is(RealMonomial)
      output = size(input)/=0
    end select
  end function
end function

! Multiplication and division by scalars.
impure elemental function multiply_RealMonomial_real(this,that) &
   & result(output)
  implicit none
  
  type(RealMonomial), intent(in) :: this
  real(dp),           intent(in) :: that
  type(RealMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient * that
end function

impure elemental function multiply_real_RealMonomial(this,that) &
   & result(output)
  implicit none
  
  real(dp),           intent(in) :: this
  type(RealMonomial), intent(in) :: that
  type(RealMonomial)             :: output
  
  output = that
  output%coefficient = output%coefficient * this
end function

impure elemental function divide_RealMonomial_real(this,that) result(output)
  implicit none
  
  type(RealMonomial), intent(in) :: this
  real(dp),           intent(in) :: that
  type(RealMonomial)             :: output
  
  output = this
  output%coefficient = output%coefficient / that
end function

! Multiplication between Monomials and Monomial-like types.
function multiply_RealMonomialable_RealMonomialable(this,that) &
   & result(output)
  implicit none
  
  class(RealMonomialable), intent(in) :: this
  class(RealMonomialable), intent(in) :: that
  type(RealMonomial)                  :: output
  
  type(RealMonomial) :: this_monomial
  type(RealMonomial) :: that_monomial
  
  integer :: i_this,i_that,i_out,ialloc
  
  this_monomial = this%to_RealMonomial()
  that_monomial = that%to_RealMonomial()
  
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
function add_RealPolynomialable_RealPolynomialable(this,that) &
   & result(output)
  implicit none
  
  class(RealPolynomialable), intent(in) :: this
  class(RealPolynomialable), intent(in) :: that
  type(RealPolynomial)                  :: output
  
  type(RealPolynomial) :: this_polynomial
  type(RealPolynomial) :: that_polynomial
  
  integer :: no_terms
  
  integer :: i,j,ialloc
  
  this_polynomial = this%to_RealPolynomial()
  that_polynomial = that%to_RealPolynomial()
  
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
    
    select type(this); type is(RealMonomial)
      select type(that); type is(RealMonomial)
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
subroutine read_RealUnivariate(this,input)
  implicit none
  
  class(RealUnivariate), intent(out) :: this
  type(String),          intent(in)  :: input
  
  type(String), allocatable :: line(:)
  integer                   :: id
  integer                   :: power
  
  select type(this); type is(RealUnivariate)
    line = split_line(input,delimiter='^')
    if (size(line)/=2) then
      call print_line(ERROR//': Unable to convert string to univariate.')
      call err()
    endif
    
    id = int(slice(line(1),2,len(line(1))))
    power = int(line(2))
    
    this = RealUnivariate(id=id,power=power)
  end select
end subroutine

function write_RealUnivariate(this) result(output)
  implicit none
  
  class(RealUnivariate), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(RealUnivariate)
    output = 'u'//this%id//'^'//this%power
  end select
end function

impure elemental function new_RealUnivariate_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(RealUnivariate)     :: this
  
  this = input
end function

subroutine read_RealMonomial(this,input)
  implicit none
  
  class(RealMonomial), intent(out) :: this
  type(String),        intent(in)  :: input
  
  type(String),         allocatable :: line(:)
  real(dp)                          :: coefficient
  type(RealUnivariate), allocatable :: modes(:)
  
  select type(this); type is(RealMonomial)
    line = split_line(input,delimiter='*')
    
    coefficient = dble(line(1))
    modes = RealUnivariate(line(2:))
    
    this = RealMonomial(coefficient, modes)
  end select
end subroutine

function write_RealMonomial(this) result(output)
  implicit none
  
  class(RealMonomial), intent(in) :: this
  type(String)                    :: output
  
  select type(this); type is(RealMonomial)
    output = this%coefficient//'*'//join(this%modes, delimiter='*')
  end select
end function

impure elemental function new_RealMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(RealMonomial)       :: this
  
  this = input
end function

subroutine read_RealPolynomial(this,input)
  implicit none
  
  class(RealPolynomial), intent(out) :: this
  type(String),          intent(in)  :: input
  
  type(String), allocatable :: terms(:)
  
  select type(this); type is(RealPolynomial)
    terms = split_line(input)
    terms = terms(filter(terms/='+'))
    this = RealPolynomial(RealMonomial(terms))
  end select
end subroutine

function write_RealPolynomial(this) result(output)
  implicit none
  
  class(RealPolynomial), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(RealPolynomial)
    output = join(this%terms, delimiter=' + ')
  end select
end function

impure elemental function new_RealPolynomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(RealPolynomial)     :: this
  
  this = input
end function
end module

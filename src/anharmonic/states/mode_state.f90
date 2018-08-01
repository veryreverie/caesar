! ======================================================================
! A generalised eigenstate along a single complex normal mode.
! ======================================================================
module mode_state_module
  use common_module
  implicit none
  
  private
  
  public :: ModeState
  public :: size
  public :: operator(+)
  public :: sum
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  
  ! An eigenstate of the form f(u)*(mw/pi)^(1/4)*e^-(w*|u|^2/2),
  !    where f is polynomial, e.g.
  !    ( a + b*(u1) + c*(u1)**2 )*(mw/pi)^(1/4)*e^(-w*|u1|^2/2)
  !    => frequency = w, coefficients = [a,b,c]
  type, extends(Stringable) :: ModeState
    integer                        :: mode_id
    real(dp)                       :: frequency
    real(dp), private, allocatable :: coefficients_(:)
  contains
    procedure, public :: coefficient => coefficient_ModeState
    procedure, public :: evaluate    => evaluate_ModeState
    
    ! I/O.
    procedure, public :: read  => read_ModeState
    procedure, public :: write => write_ModeState
  end type
  
  interface ModeState
    module procedure new_ModeState
    module procedure new_ModeState_mode
    module procedure new_ModeState_String
  end interface
  
  interface size
    module procedure size_ModeState
  end interface
  
  interface operator(+)
    module procedure add_ModeState_ModeState
  end interface
  
  interface sum
    module procedure sum_ModeStates
  end interface
  
  interface operator(-)
    module procedure subtract_ModeState_ModeState
  end interface
  
  interface operator(*)
    module procedure multiply_ModeState_real
    module procedure multiply_real_ModeState
  end interface
  
  interface operator(/)
    module procedure divide_ModeState_real
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality:
!    - constructor
!    - size() function
!    - coefficient() getter
! ----------------------------------------------------------------------
function new_ModeState(mode_id,frequency,coefficients) result(this)
  implicit none
  
  integer,  intent(in) :: mode_id
  real(dp), intent(in) :: coefficients(:)
  real(dp), intent(in) :: frequency
  type(ModeState)      :: this
  
  this%mode_id       = mode_id
  this%frequency     = frequency
  this%coefficients_ = coefficients
end function

function new_ModeState_mode(mode,frequency,coefficients) result(this)
  implicit none
  
  type(ComplexMode), intent(in)           :: mode
  real(dp),          intent(in), optional :: frequency
  real(dp),          intent(in)           :: coefficients(:)
  type(ModeState)                         :: this
  
  this%mode_id = mode%id
  if (present(frequency)) then
    this%frequency = frequency
  else
    this%frequency = mode%frequency
  endif
  this%coefficients_ = coefficients
end function

! Returns the largest power with a non-zero coefficient.
function size_ModeState(this) result(output)
  implicit none
  
  type(ModeState), intent(in) :: this
  integer                     :: output
  
  output = size(this%coefficients_)-1
end function

! coefficient(i) returns the coefficient of u^i.
! Intended to avoid all the +/-1s when indexing,
!    since this%coefficients_(i) is the coefficient of u^(i-1)
function coefficient_ModeState(this,power) result(output)
  implicit none
  
  class(ModeState), intent(in) :: this
  integer,          intent(in) :: power
  real(dp)                     :: output
  
  output = this%coefficients_(power+1)
end function

! ----------------------------------------------------------------------
! Arithmetic operations on states.
! ----------------------------------------------------------------------
function add_ModeState_ModeState(this,that) &
   & result(output)
  implicit none
  
  type(ModeState), intent(in) :: this
  type(ModeState), intent(in) :: that
  type(ModeState)             :: output
  
  integer :: max_coefficient
  
  integer :: i,ialloc
  
  ! Check inputs are consistent.
  if (this%mode_id/=that%mode_id) then
    call print_line(CODE_ERROR//': Trying to add together states along &
       &diferent modes.')
    call err()
  endif
  
  ! Copy across everything but coefficients from one input.
  output = this
  
  ! Add together the coefficients of both inputs.
  max_coefficient = max(size(this),size(that))
  allocate( output%coefficients_(max_coefficient+1), &
          & stat=ialloc); call err(ialloc)
  output%coefficients_ = 0
  do i=0,size(this)
    output%coefficients_(i+1) = output%coefficients_(i+1) &
                            & + this%coefficient(i)
  enddo
  
  do i=0,size(that)
    output%coefficients_(i+1) = output%coefficients_(i+1) &
                            & + that%coefficient(i)
  enddo
end function

function sum_ModeStates(input) result(output)
  implicit none
  
  type(ModeState), intent(in) :: input(:)
  type(ModeState)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(ERROR//': Trying to sum an empty array.')
    call err()
  endif
  
  output = input(1)
  
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

function subtract_ModeState_ModeState(this,that) &
   & result(output)
  implicit none
  
  type(ModeState), intent(in) :: this
  type(ModeState), intent(in) :: that
  type(ModeState)             :: output
  
  integer :: max_coefficient
  
  integer :: i,ialloc
  
  ! Check inputs are consistent.
  if (this%mode_id/=that%mode_id) then
    call print_line(CODE_ERROR//': Trying to add together states along &
       &diferent modes.')
    call err()
  endif
  
  ! Copy across everything but coefficients from one input.
  output = this
  
  ! Subtract the coefficients of that from those of this.
  max_coefficient = max(size(this),size(that))
  allocate( output%coefficients_(max_coefficient+1), &
          & stat=ialloc); call err(ialloc)
  output%coefficients_ = 0
  do i=0,size(this)
    output%coefficients_(i+1) = output%coefficients_(i+1) &
                            & + this%coefficient(i)
  enddo
  
  do i=0,size(that)
    output%coefficients_(i+1) = output%coefficients_(i+1) &
                            & - that%coefficient(i)
  enddo
end function

function multiply_ModeState_real(this,that) result(output)
  implicit none
  
  type(ModeState), intent(in) :: this
  real(dp),        intent(in) :: that
  type(ModeState)             :: output
  
  output = this
  output%coefficients_ = output%coefficients_ * that
end function

function multiply_real_ModeState(this,that) result(output)
  implicit none
  
  real(dp),        intent(in) :: this
  type(ModeState), intent(in) :: that
  type(ModeState)             :: output
  
  output = that
  output%coefficients_ = output%coefficients_ * this
end function

function divide_ModeState_real(this,that) result(output)
  implicit none
  
  type(ModeState), intent(in) :: this
  real(dp),        intent(in) :: that
  type(ModeState)             :: output
  
  output = this
  output%coefficients_ = output%coefficients_ / that
end function

! ----------------------------------------------------------------------
! Evaluates the state at a given displacement, u, along the normal mode.
! ----------------------------------------------------------------------
! N.B. does not include the normalisation factor of (m*freq/pi)^(1/4).
function evaluate_ModeState(this,u) result(output)
  implicit none
  
  class(ModeState),                intent(in) :: this
  type(ComplexSingleDisplacement), intent(in) :: u
  real(dp)                                    :: output
  
  real(dp) :: term
  
  integer :: i
  
  output = 0
  term = 1
  do i=0,size(this)
    output = output + this%coefficient(i) * term
    term = term*u%magnitude
  enddo
  output = output * exp(-0.5_dp*this%frequency*abs(u%magnitude*u%magnitude))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
! A state with coefficients=[3.1e2,7.4e-2,-4.9e-6] and
!    frequency=[6.6e-4] along mode 7 becomes
! 'u7: (3.1e2 +7.4e-2*u^1 -4.9e-6*u^2)*(mw/pi)^(1/4)*e^(-6.6e-4*|u|^2)'
subroutine read_ModeState(this,input)
  implicit none
  
  class(ModeState), intent(out) :: this
  type(String),     intent(in)  :: input
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: term(:)
  
  integer               :: mode_id
  real(dp)              :: frequency
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ModeState)
    ! Split off the mode id term.
    line = split_line(input,delimiter=':')
    mode_id = int(slice(line(1),2,len(line(1))))
    
    ! Split into coefficients and exponent.
    ! line = [ '(3.1e2 +7.4e-2*u^1 -4.9e-6*u^2',
    !          '*(mw/pi',
    !          '^(1/4
    !          '*e^(-6.6e-4*|u|^2'                ]
    line = split_line(line(2),delimiter=')')
    if (size(line)/=4) then
      call print_line(ERROR//': Unable to parse string into ModeState:')
      call print_line(input)
      call err()
    endif
    
    ! Extract the frequency from 'e^(-6.6e-4*|u|^2'.
    ! Trims the '*e^(' from the front, and the '*|u|^2' from the back,
    !    then finds the negative of the remaining real.
    frequency = -dble(slice(line(4),5,len(line(2))-6))
    
    ! Extract the coefficients from '(3.1e2 +7.4e-2*u^1 -4.9e-6*u^2'.
    ! Trim the '(', and split the string by spaces.
    ! line = ['3.1e2', '+7.4e-2*u^1', '-4.9e-6*u^2']
    line = split_line(slice(line(1),2,len(line(1))))
    allocate(coefficients(size(line)),stat=ialloc); call err(ialloc)
    do i=1,size(coefficients)
      term = split_line(line(i),delimiter='*')
      coefficients(i) = dble(term(1))
    enddo
    
    this = ModeState(mode_id, frequency, coefficients)
  class default
    call err()
  end select
end subroutine

function write_ModeState(this) result(output)
  implicit none
  
  class(ModeState), intent(in) :: this
  type(String)                 :: output
  
  integer :: i
  
  select type(this); type is(ModeState)
    output = 'u'//this%mode_id//': ('//this%coefficient(0)
    do i=1,size(this)
      if (this%coefficient(i)>=0) then
        output = output//' +'//abs(this%coefficient(i))
      else
        output = output//' -'//abs(this%coefficient(i))
      endif
      output = output//'*u^'//i
    enddo
    output = output//')*(mw/pi)^(1/4)*e^('
    if (this%frequency>=0) then
      output = output//'-'//abs(this%frequency)//'*|u|^2/2)'
    else
      output = output//'+'//abs(this%frequency)//'*|u|^2/2)'
    endif
  class default
    call err()
  end select
end function

impure elemental function new_ModeState_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ModeState)          :: this
  
  call this%read(input)
end function
end module

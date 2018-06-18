! ======================================================================
! A generalised eigenstate along a single complex normal mode.
! ======================================================================
module single_mode_state_module
  use common_module
  implicit none
  
  private
  
  public :: SingleModeState
  public :: size
  public :: operator(+)
  public :: sum
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: generate_harmonic_basis
  
  ! An eigenstate of the form f(u)*e^-(w*u^2/2) where f is polynomial, e.g.
  !    ( a + b*(u1) + c*(u1)**2 )*e^(-w*(u1)^2/2) => frequency    = w,
  !                                                  coefficients = [a,b,c]
  type, extends(Stringable) :: SingleModeState
    integer                        :: mode_id
    real(dp)                       :: frequency
    real(dp), private, allocatable :: coefficients_(:)
  contains
    procedure, public :: coefficient => coefficient_SingleModeState
    procedure, public :: evaluate    => evaluate_SingleModeState
    procedure, public :: bloch_momentum
    
    ! I/O.
    procedure, public :: read  => read_SingleModeState
    procedure, public :: write => write_SingleModeState
  end type
  
  interface SingleModeState
    module procedure new_SingleModeState
    module procedure new_SingleModeState_String
  end interface
  
  interface size
    module procedure size_SingleModeState
  end interface
  
  interface operator(+)
    module procedure add_SingleModeState_SingleModeState
  end interface
  
  interface sum
    module procedure sum_SingleModeStates
  end interface
  
  interface operator(-)
    module procedure subtract_SingleModeState_SingleModeState
  end interface
  
  interface operator(*)
    module procedure multiply_SingleModeState_real
    module procedure multiply_real_SingleModeState
  end interface
  
  interface operator(/)
    module procedure divide_SingleModeState_real
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality:
!    - constructor
!    - size() function
!    - coefficient() getter
! ----------------------------------------------------------------------
function new_SingleModeState(complex_mode,coefficients,frequency) result(this)
  implicit none
  
  type(ComplexMode), intent(in)           :: complex_mode
  real(dp),          intent(in)           :: coefficients(:)
  real(dp),          intent(in), optional :: frequency
  type(SingleModeState)                   :: this
  
  this%mode_id = complex_mode%id
  if (present(frequency)) then
    this%frequency = frequency
  else
    this%frequency = complex_mode%frequency
  endif
  this%coefficients_ = coefficients
end function

! Returns the largest power with a non-zero coefficient.
function size_SingleModeState(this) result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: this
  integer                           :: output
  
  output = size(this%coefficients_)-1
end function

! coefficient(i) returns the coefficient of u^i.
! Intended to avoid all the +/-1s when indexing,
!    since this%coefficients_(i) is the coefficient of u^(i-1)
function coefficient_SingleModeState(this,power) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  integer,                intent(in) :: power
  real(dp)                           :: output
  
  output = this%coefficients_(power+1)
end function

! ----------------------------------------------------------------------
! Arithmetic operations on states.
! ----------------------------------------------------------------------
function add_SingleModeState_SingleModeState(this,that) &
   & result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: this
  type(SingleModeState), intent(in) :: that
  type(SingleModeState)             :: output
  
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

function sum_SingleModeStates(input) result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: input(:)
  type(SingleModeState)             :: output
  
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

function subtract_SingleModeState_SingleModeState(this,that) &
   & result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: this
  type(SingleModeState), intent(in) :: that
  type(SingleModeState)             :: output
  
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

function multiply_SingleModeState_real(this,that) result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: this
  real(dp),              intent(in) :: that
  type(SingleModeState)             :: output
  
  output = this
  output%coefficients_ = output%coefficients_ * that
end function

function multiply_real_SingleModeState(this,that) result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(SingleModeState), intent(in) :: that
  type(SingleModeState)             :: output
  
  output = that
  output%coefficients_ = output%coefficients_ * this
end function

function divide_SingleModeState_real(this,that) result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: this
  real(dp),              intent(in) :: that
  type(SingleModeState)             :: output
  
  output = this
  output%coefficients_ = output%coefficients_ / that
end function

! ----------------------------------------------------------------------
! Generates the harmonic basis functions along a specific mode.
! ----------------------------------------------------------------------
! Uses the recurrence relation:
!   |n> = sqrt(2*freq/n) u |n-1> - sqrt((n-1)/n) |n-2>
! N.B. basis_functions(i) = |i-1> because |0> is a state.
function generate_harmonic_basis(mode,harmonic_states_cutoff,frequency) &
   & result(output)
  implicit none
  
  type(ComplexMode), intent(in)           :: mode
  integer,           intent(in)           :: harmonic_states_cutoff
  real(dp),          intent(in), optional :: frequency
  type(SingleModeState), allocatable      :: output(:)
  
  real(dp)              :: freq
  real(dp)              :: normalisation
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,ialloc
  
  if (present(frequency)) then
    freq = frequency
  else
    freq = mode%frequency
  endif
  
  normalisation = (freq/PI)**0.25_dp
  
  allocate(output(harmonic_states_cutoff+1), stat=ialloc); call err(ialloc)
  
  ! Calculate |0>.
  if (harmonic_states_cutoff >= 0) then
    coefficients = [normalisation]
    output(1) = SingleModeState(mode,coefficients,frequency)
  endif
  
  ! Calculate |1>.
  if (harmonic_states_cutoff >= 1) then
    coefficients = [0.0_dp, normalisation*sqrt(freq*2)]
    output(2) = SingleModeState(mode,coefficients,frequency)
  endif
  
  do i=3,size(output)
    ! Calculate |i-1> from |i-2> and |i-3>.
    deallocate(coefficients, stat=ialloc); call err(ialloc)
    allocate(coefficients(i), stat=ialloc); call err(ialloc)
    coefficients = 0
    
    ! sqrt(2*freq/(i-1)) u |i-2>.
    coefficients(2:) = coefficients(2:)        &
                   & + sqrt(2*frequency/(i-1)) &
                   & * output(i-1)%coefficients_(:i-1)
    
    ! -sqrt((i-2)/(i-1)) |i-3>.
    coefficients(:i-2) = coefficients(:i-2)     &
                     & - sqrt((i-2.0_dp)/(i-1)) &
                     & * output(i-2)%coefficients_
    
    output(i) = SingleModeState(mode,coefficients,frequency)
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates the state at a given displacement, u, along the normal mode.
! ----------------------------------------------------------------------
function evaluate_SingleModeState(this,u) result(output)
  implicit none
  
  class(SingleModeState),        intent(in) :: this
  type(ComplexSingleModeVector), intent(in) :: u
  real(dp)                                  :: output
  
  real(dp) :: term
  
  integer :: i
  
  output = 0
  term = 1
  do i=0,size(this)
    output = output + this%coefficient(i) * term
    term = term*u%magnitude
  enddo
  output = output * exp(-0.5_dp*this%frequency*u%magnitude*u%magnitude)
end function

! ----------------------------------------------------------------------
! Returns the Bloch momentum of the state.
! ----------------------------------------------------------------------
! The state |i> at q-point q has a Bloch momentum of i*q.
function bloch_momentum(this,mode,qpoints) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  type(ComplexMode),      intent(in) :: mode
  type(QpointData),       intent(in) :: qpoints(:)
  type(FractionVector)               :: output
  
  type(QpointData) :: qpoint
  
  if (this%mode_id/=mode%id) then
    call print_line(CODE_ERROR//': Trying to find the Bloch momentum of an &
       &incompatible state and mode.')
    call err()
  endif
  
  qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
  
  output = qpoint%qpoint * size(this)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
! A state with coefficients=[3.1e2,7.4e-2,-4.9e-6] and
!    frequency=[6.6e-4] becomes
! '(3.1e2 +7.4e-2*u^1 -4.9e-6*u^2)e^(-6.6e-4*|u|^2)'
subroutine read_SingleModeState(this,input)
  implicit none
  
  class(SingleModeState), intent(out) :: this
  type(String),           intent(in)  :: input
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: term(:)
  
  real(dp)              :: frequency
  real(dp), allocatable :: coefficients(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SingleModeState)
    ! Split the input into coefficients and exponent.
    ! line = ['(3.1e2 +7.4e-2*u^1 -4.9e-6*u^2', 'e^(-6.6e-4*|u|^2']
    line = split_line(input,delimiter=')')
    if (size(line)/=2) then
      call print_line(ERROR//': Unable to parse string into SingleModeState:')
      call print_line(input)
      call err()
    endif
    
    ! Extract the frequency from 'e^(-6.6e-4*|u|^2'.
    ! Trims the 'e^(' from the front, and the '*|u|^2' from the back,
    !    then finds the negative of the remaining real.
    frequency = -dble(slice(line(2),4,len(line(2))-8))
    
    ! Split the coefficients string by spaces.
    ! line = ['3.1e2', '+7.4e-2*u^1', '-4.9e-6*u^2']
    line = split_line(slice(line(1),2,len(line(1))))
    allocate(coefficients(size(line)),stat=ialloc); call err(ialloc)
    do i=1,size(coefficients)
      term = split_line(line(i),delimiter='*')
      coefficients(i) = dble(term(1))
    enddo
    
    this = SingleModeState(frequency, coefficients)
  end select
end subroutine

function write_SingleModeState(this) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  type(String)                       :: output
  
  integer :: i
  
  select type(this); type is(SingleModeState)
    output = '('//this%coefficient(0)
    do i=1,size(this)
      if (this%coefficient(i)>=0) then
        output = output//' +'//abs(this%coefficient(i))
      else
        output = output//' -'//abs(this%coefficient(i))
      endif
      output = output//'*u^'//i
    enddo
    output = output//')e^('
    if (this%frequency>=0) then
      output = output//'-'//abs(this%frequency)//'*|u|^2/2)'
    else
      output = output//'+'//abs(this%frequency)//'*|u|^2/2)'
    endif
  end select
end function

impure elemental function new_SingleModeState_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(SingleModeState)    :: this
  
  this = input
end function
end module

! ======================================================================
! A generalised eigenstate along a single normal mode.
! ======================================================================
module single_mode_states_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! An eigenstate of the form f(u)*e^-(w*u^2/2) where f is polynomial, e.g.
  !    ( a + b*(u1) + c*(u1)**2 )*e^(-w*(u1)^2/2) => frequency    = w,
  !                                                  coefficients = [a,b,c]
  type, extends(Stringable) :: SingleModeState
    integer                        :: mode
    real(dp)                       :: frequency
    real(dp), private, allocatable :: coefficients_(:)
  contains
    procedure, public :: coefficient => coefficient_SingleModeState
    
    generic, public :: operator(+) => add_SingleModeState_SingleModeState
    generic, public :: operator(-) => subtract_SingleModeState_SingleModeState
    generic, public :: operator(*) => multiply_SingleModeState_real, &
                                    & multiply_real_SingleModeState
    generic, public :: operator(/) => divide_SingleModeState_real
    
    procedure :: evaluate => evaluate_SingleModeState
    
    procedure, private             :: add_SingleModeState_SingleModeState
    procedure, private             :: subtract_SingleModeState_SingleModeState
    procedure, private             :: multiply_SingleModeState_real
    procedure, private, pass(that) :: multiply_real_SingleModeState
    procedure, private             :: divide_SingleModeState_real
    
    ! I/O.
    procedure :: str => str_SingleModeState
  end type
  
  interface size
    module procedure size_SingleModeState
  end interface
contains

! ----------------------------------------------------------------------
! Returns the largest power of a single mode state.
! ----------------------------------------------------------------------
function size_SingleModeState(this) result(output)
  implicit none
  
  type(SingleModeState), intent(in) :: this
  integer                           :: output
  
  output = size(this%coefficients_)-1
end function

! ----------------------------------------------------------------------
! coefficient(i) returns the coefficient of u^i.
! ----------------------------------------------------------------------
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
  
  class(SingleModeState), intent(in) :: this
  class(SingleModeState), intent(in) :: that
  type(SingleModeState)              :: output
  
  integer :: i,ialloc
  
  if (this%mode/=that%mode) then
    call err()
  endif
  
  output%mode = this%mode
  output%frequency = this%frequency
  
  allocate( output%coefficients_(max(size(this),size(that))+1), &
          & stat=ialloc); call err(ialloc)
  output%coefficients_ = 0
  
  do i=0,size(this)
    output%coefficients_(i+1) = output%coefficients_(i+1) + this%coefficient(i)
  enddo
  
  do i=0,size(that)
    output%coefficients_(i+1) = output%coefficients_(i+1) + that%coefficient(i)
  enddo
end function

function subtract_SingleModeState_SingleModeState(this,that) &
   & result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  class(SingleModeState), intent(in) :: that
  type(SingleModeState)              :: output
  
  integer :: i,ialloc
  
  if (this%mode/=that%mode) then
    call err()
  endif
  
  output%mode = this%mode
  output%frequency = this%frequency
  
  allocate( output%coefficients_(max(size(this),size(that))+1), &
          & stat=ialloc); call err(ialloc)
  output%coefficients_ = 0
  
  do i=0,size(this)
    output%coefficients_(i+1) = output%coefficients_(i+1) + this%coefficient(i)
  enddo
  
  do i=0,size(that)
    output%coefficients_(i+1) = output%coefficients_(i+1) - that%coefficient(i)
  enddo
end function

function multiply_SingleModeState_real(this,that) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  real(dp),               intent(in) :: that
  type(SingleModeState)              :: output
  
  output = this
  output%coefficients_ = output%coefficients_ * that
end function

function multiply_real_SingleModeState(this,that) result(output)
  implicit none
  
  real(dp),               intent(in) :: this
  class(SingleModeState), intent(in) :: that
  type(SingleModeState)              :: output
  
  output = that
  output%coefficients_ = output%coefficients_ * this
end function

function divide_SingleModeState_real(this,that) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  real(dp),               intent(in) :: that
  type(SingleModeState)              :: output
  
  output = this
  output%coefficients_ = output%coefficients_ / that
end function

! ----------------------------------------------------------------------
! I/O on states.
! ----------------------------------------------------------------------
function str_SingleModeState(this) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  type(String)                       :: output
  
  type(String) :: token
  integer      :: i
  
  output = '('
  do i=0,size(this)
    token = trim(str(this%coefficient(i)))
    if (slice(token,1,1)=='-' .or. i==0) then
      output = output//token
    else
      output = output//'+'//token
    endif
    output = output//'u^'//i
  enddo
  
  output = output//')e^('
  token = trim(str(this%frequency))
  if (slice(token,1,1)=='-') then
    output = output//'+'
    token = slice(token,2,len(token))
  else
    output = output//'-'
  endif
  output = output//token//'*u^2)'
end function

! ----------------------------------------------------------------------
! Evaluates the state at a given displacement, u, along the normal mode.
! ----------------------------------------------------------------------
function evaluate_SingleModeState(this,u) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  real(dp),               intent(in) :: u
  real(dp)                           :: output
  
  integer :: i
  
  output = 0
  do i=0,size(this)
    output = output + this%coefficient(i) * u**i
  enddo
  output = output * exp(-0.5_dp*this%frequency*u*u)
end function

! ----------------------------------------------------------------------
! Generates the harmonic basis functions along a specific mode.
! ----------------------------------------------------------------------
! Uses the recurrence relation:
!   |n> = sqrt(2*freq/n) u |n-1> - sqrt((n-1)/n) |n-2>
! N.B. basis_functions(i) = |i-1> because |0> is a state.
function generate_harmonic_basis(mode,frequency,harmonic_states_cutoff) &
   & result(output)
  use constants_module, only : pi
  use coupling_module
  implicit none
  
  integer,  intent(in)               :: mode
  real(dp), intent(in)               :: frequency
  integer,  intent(in)               :: harmonic_states_cutoff
  type(SingleModeState), allocatable :: output(:)
  
  real(dp) :: normalisation
  
  integer :: i,ialloc
  
  normalisation = (frequency/pi)**0.25_dp
  
  allocate(output(harmonic_states_cutoff+1), stat=ialloc); call err(ialloc)
  
  ! Every basis function has the same exponent.
  do i=1,size(output)
    output(i)%mode = mode
    output(i)%frequency = frequency
  enddo
  
  ! Calculate |0>.
  if (harmonic_states_cutoff >= 0) then
    output(1)%coefficients_ = [normalisation]
  endif
  
  ! Calculate |1>.
  if (harmonic_states_cutoff >= 1) then
    output(2)%coefficients_ = [0.0_dp, normalisation*sqrt(frequency*2)]
  endif
  
  do i=3,size(output)
    ! Calculate |i-1> from |i-2> and |i-3>.
    allocate(output(i)%coefficients_(i), stat=ialloc); call err(ialloc)
    output(i)%coefficients_ = 0
    ! sqrt(2*freq/(i-1)) u |i-2>.
    output(i)%coefficients_(2:) = output(i)%coefficients_(2:) &
                              & + sqrt(2*frequency/(i-1))     &
                              & * output(i-1)%coefficients_(:i-1)
    ! -sqrt((i-2)/(i-1)) |i-3>.
    output(i)%coefficients_(:i-2) = output(i)%coefficients_(:i-2) &
                                & - sqrt((i-2.0_dp)/(i-1))        &
                                & * output(i-2)%coefficients_
  enddo
end function
end module

! ======================================================================
! A generalised eigenstate basis function.
! ======================================================================
module eigenstates_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! An eigenstate of the form f(u)*e^-(w*u^2/2) where f is polynomial, e.g.
  !    ( a + b*(u1) + c*(u1)**2 )*e^(-w*(u1)^2/2) => frequency    = w,
  !                                                  coefficients = [a,b,c]
  type :: SingleModeState
    real(dp)              :: frequency
    real(dp), allocatable :: coefficients(:)
  contains
    procedure :: evaluate => evaluate_SingleModeState_real
  end type
contains

! ------------------------------------------------------------
! Evaluates the basis function at a given displacement along the normal mode.
! ------------------------------------------------------------
function evaluate_SingleModeState_real(this,u) result(output)
  implicit none
  
  class(SingleModeState), intent(in) :: this
  real(dp),               intent(in) :: u
  real(dp)                           :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this%coefficients)
    output = output + this%coefficients(i)*u**(i-1)
  enddo
  output = output * exp(-0.5_dp*this%frequency*u*u)
end function

! ------------------------------------------------------------
! Generates the harmonic basis functions along a specific mode.
! ------------------------------------------------------------
! Uses the recurrence relation:
!   |n> = sqrt(2*freq/n) |n-1> - sqrt((n-1)/n) |n-2>
! N.B. basis_functions(i) = |i-1> because |0> is a state.
function generate_harmonic_basis(frequency,no_harmonic_states) &
   & result(output)
  use constants_module, only : pi
  use coupling_module
  implicit none
  
  real(dp), intent(in)               :: frequency
  integer,  intent(in)               :: no_harmonic_states
  type(SingleModeState), allocatable :: output(:)
  
  real(dp) :: normalisation
  
  integer :: i,ialloc
  
  normalisation = (frequency/pi)**0.25_dp
  
  allocate(output(no_harmonic_states), stat=ialloc); call err(ialloc)
  
  ! Every basis function has the same exponent.
  do i=1,no_harmonic_states
    output(i)%frequency = frequency
  enddo
  
  ! Calculate |0>.
  if (no_harmonic_states >= 1) then
    output(1)%coefficients = [normalisation]
  endif
  
  ! Calculate |1>.
  if (no_harmonic_states >= 2) then
    output(2)%coefficients = [0.0_dp, normalisation*sqrt(frequency*2)]
  endif
  
  do i=3,no_harmonic_states
    ! Calculate |i-1> from |i-2> and |i-3>.
    allocate(output(i)%coefficients(i), stat=ialloc); call err(ialloc)
    output(i)%coefficients = 0
    ! sqrt(2*freq/(i-1)) u |i-2>.
    output(i)%coefficients(2:) = output(i)%coefficients(2:) &
                             & + sqrt(2*frequency/(i-1))    &
                             & * output(i-1)%coefficients(:i-1)
    ! -sqrt((i-2)/(i-1)) |i-3>.
    output(i)%coefficients(:i-2) = output(i)%coefficients(:i-2) &
                               & - sqrt((i-2.0_dp)/(i-1))       &
                               & * output(i-2)%coefficients
  enddo
end function
end module

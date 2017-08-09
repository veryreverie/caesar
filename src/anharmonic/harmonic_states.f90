! ======================================================================
! A generalised harmonic basis function.
! ======================================================================
module harmonic_states_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! A harmonic eigenstate, e.g.
  !    ( a + b(u1) + c(u1)**2 )*e**(E*(u1)**2) => exponent=E, coeffs=[a,b,c]
  type :: HarmonicState
    real(dp)              :: frequency
    real(dp), allocatable :: coefficients(:)
  contains
    procedure :: evaluate => evaluate_HarmonicState_real
  end type
contains

! ------------------------------------------------------------
! Evaluates the basis function at a given displacement along the normal mode.
! ------------------------------------------------------------
function evaluate_HarmonicState_real(this,u) result(output)
  implicit none
  
  class(HarmonicState), intent(in) :: this
  real(dp),             intent(in) :: u
  real(dp)                         :: output
  
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
  
  real(dp), intent(in)             :: frequency
  integer,  intent(in)             :: no_harmonic_states
  type(HarmonicState), allocatable :: output(:)
  
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

! ----------------------------------------------------------------------
! Calculates integrals of the form I(n) = u^n*e^(-freq*u*u)
! ----------------------------------------------------------------------
! N.B. output(i) = I(i-1) because I(0) exists.
! Uses the relations:
!    I(n) = integral[ u^n*e^(-freq*u*u) ] from -infinity to infinity.
!    I(0) = sqrt(pi/freq)
!    I(1) = 0
!    I(n) = (n-1)/(2*freq) * I(n-2)
function calculate_gaussian_integrals(frequency,no_harmonic_states, &
   & no_basis_functions) result(output)
  use constants_module, only : pi
  implicit none
  
  real(dp), intent(in)  :: frequency
  integer,  intent(in)  :: no_harmonic_states
  integer,  intent(in)  :: no_basis_functions
  real(dp), allocatable :: output(:)
  
  integer :: i,ialloc
  
  allocate( output(2*no_harmonic_states+no_basis_functions), &
          & stat=ialloc); call err(ialloc)
  output(1) = sqrt(pi/frequency)
  output(2) = 0
  do i=3,size(output)
    output(i) = output(i-2)*i/(2*frequency)
  enddo
end function

! ----------------------------------------------------------------------
! Generates all elements <i|u^k|j> along a given mode.
! ----------------------------------------------------------------------
! output(k) is a matrix whose i,j element is <i-1|u^(k-1)|j-1>.
! The -1 offsets are due to |0> and u^0.
function generate_harmonic_couplings(harmonic_states,no_basis_functions, &
   & gaussian_integrals) result(output)
  use constants_module, only : pi
  use linear_algebra_module
  implicit none
  
  type(HarmonicState), intent(in) :: harmonic_states(:)
  integer,             intent(in) :: no_basis_functions
  real(dp),            intent(in) :: gaussian_integrals(:)
  type(RealMatrix), allocatable   :: output(:)
  
  real(dp), allocatable :: matrix(:,:)
  
  integer :: i,j,k,i2,j2,ialloc
  
  allocate( matrix(size(harmonic_states),size(harmonic_states)), &
          & output(no_basis_functions),                          &
          & stat=ialloc); call err(ialloc)
  
  do k=1,no_basis_functions
    ! Calculate output(k) = {<i|u^k|j>}
    matrix = 0
    do j=1,size(harmonic_states)
      do i=1,j
        do j2=1,size(harmonic_states(j)%coefficients)
          do i2=1,size(harmonic_states(i)%coefficients)
            ! N.B. the -2 is because i2 refers to u^(i2-1) etc.
            matrix(i,j) = matrix(i,j)                         &
                      & + harmonic_states(j)%coefficients(j2) &
                      & * harmonic_states(i)%coefficients(i2) &
                      & * gaussian_integrals(i2+j2+k-2)
          enddo
        enddo
        matrix(j,i) = matrix(i,j)
      enddo
    enddo
    output(k) = matrix
  enddo
end function
end module

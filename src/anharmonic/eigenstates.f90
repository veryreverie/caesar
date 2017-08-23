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
    procedure :: evaluate => evaluate_SingleModeState
  end type
  
  ! A product of single mode states.
  type :: ProductState
    type(SingleModeState), allocatable :: states(:)
  contains
    procedure :: evaluate => evaluate_ProductState
  end type
contains

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
  do i=1,size(this%coefficients)
    output = output + this%coefficients(i)*u**(i-1)
  enddo
  output = output * exp(-0.5_dp*this%frequency*u*u)
end function

! ----------------------------------------------------------------------
! Evaluates the state at a given displacement in normal mode co-ordinates.
! ----------------------------------------------------------------------
function evaluate_ProductState(this,displacement) result(output)
  use normal_mode_module
  implicit none
  
  class(ProductState), intent(in) :: this
  type(ModeVector),    intent(in) :: displacement
  real(dp)                        :: output
  
  integer :: i
  
  output = 1
  do i=1,size(this%states)
    output = output*this%states(i)%evaluate(displacement%vector(i))
  enddo
end function

! ----------------------------------------------------------------------
! Generates the harmonic basis functions along a specific mode.
! ----------------------------------------------------------------------
! Uses the recurrence relation:
!   |n> = sqrt(2*freq/n) |n-1> - sqrt((n-1)/n) |n-2>
! N.B. basis_functions(i) = |i-1> because |0> is a state.
function generate_harmonic_basis(frequency,harmonic_states_cutoff) &
   & result(output)
  use constants_module, only : pi
  use coupling_module
  implicit none
  
  real(dp), intent(in)               :: frequency
  integer,  intent(in)               :: harmonic_states_cutoff
  type(SingleModeState), allocatable :: output(:)
  
  real(dp) :: normalisation
  
  integer :: i,ialloc
  
  normalisation = (frequency/pi)**0.25_dp
  
  allocate(output(harmonic_states_cutoff+1), stat=ialloc); call err(ialloc)
  
  ! Every basis function has the same exponent.
  do i=1,harmonic_states_cutoff
    output(i)%frequency = frequency
  enddo
  
  ! Calculate |0>.
  if (harmonic_states_cutoff >= 1) then
    output(1)%coefficients = [normalisation]
  endif
  
  ! Calculate |1>.
  if (harmonic_states_cutoff >= 2) then
    output(2)%coefficients = [0.0_dp, normalisation*sqrt(frequency*2)]
  endif
  
  do i=3,size(output)
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
! Constructs all relevant product states from single-mode states.
! ----------------------------------------------------------------------
function construct_product_states(single_mode_states,coupling) result(output)
  use coupling_module
  use grid_types_module
  implicit none
  
  type(SingleModeState), intent(in) :: single_mode_states(:,:)
  type(CoupledModes),    intent(in) :: coupling(:)
  type(ProductState), allocatable   :: output(:)
  
  integer              :: harmonic_states_cutoff
  integer              :: no_modes
  integer              :: output_size
  integer, allocatable :: grid(:,:)
  
  integer :: i,j,k,l,ialloc
  
  ! States range from |0> to |n>, so n=size(states)-1.
  harmonic_states_cutoff = size(single_mode_states,1)-1
  no_modes = size(single_mode_states,2)
  
  ! Calculate how much space is needed and allocate output.
  output_size = 0
  do i=1,size(coupling)
    ! This gives the array length of 'grid' below.
    output_size = output_size &
              & + octahedral_grid_size( size(coupling(i)),        &
              &                         harmonic_states_cutoff-1, &
              &                         include_negatives=.false.)
  enddo
  allocate(output(output_size), stat=ialloc); call err(ialloc)
  
  ! Calculate output, coupling by coupling.
  j = 0
  do i=1,size(coupling)
    ! Generate an octahedral grid of points.
    ! Each poing will be mapped on to the states in each product.
    ! The -1 comes from not taking any |0> states, since they will be handled
    !    by subsidiary couplings.
    grid = generate_octahedral_grid( size(coupling(i)),        &
                                   & harmonic_states_cutoff-1, &
                                   & include_negatives=.false.)
    do k=1,size(grid)
      allocate(output(j+k)%states(no_modes), stat=ialloc); call err(ialloc)
      ! Initialise all single-mode states to |0>.
      do l=1,no_modes
        output(j+k)%states(l) = single_mode_states(1,l)
      enddo
      ! Set all states in the coupling.
      ! One +1 in the +2 comes for the same reason as the -1 above.
      ! The other +1 comes from single_mode_states(1,i) = |0>.
      do l=1,size(coupling(i))
        output(j+k)%states(coupling(i)%modes(l)) = &
           & single_mode_states(grid(l,k)+2,coupling(i)%modes(l))
      enddo
    enddo
    j = j+size(grid)
  enddo
end function
end module

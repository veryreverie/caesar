! ======================================================================
! Calculates corrections due to perturbation theory.
! ======================================================================
module perturbation_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  private
  
  ! The energy and state corrections when a perturbation is real.
  type, public :: RealPerturbation
    real(dp), allocatable :: energy(:)
    real(dp), allocatable :: states(:,:)
  end type
  
  ! The energy and state corrections when a perturbation is complex.
  type, public :: ComplexPerturbation
    real(dp), allocatable :: energy(:)
    complex(dp), allocatable :: states(:,:)
  end type
  
  public :: calculate_perturbation
  
  interface calculate_perturbation
    module procedure calculate_perturbation_real
    module procedure calculate_perturbation_complex
  end interface
contains

! Calculates corrections to the energy.
function calculate_energy_correction(eigenvalues,perturbation, &
   & energy_corrections,state_corrections) result(output)
  implicit none
  
  real(dp), intent(in)           :: eigenvalues(:)
  real(dp), intent(in)           :: perturbation(:,:)
  real(dp), intent(in), optional :: energy_corrections(:,:)
  real(dp), intent(in), optional :: state_corrections(:,:,:)
  real(dp), allocatable          :: output(:)
  
  ! Previous corrections.
  real(dp),    allocatable :: energy(:,:)
  real(dp), allocatable :: states(:,:,:)
  
  ! Array lengths.
  integer :: no_states
  integer :: order
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  no_states = size(eigenvalues)
  
  if (present(energy_corrections)) then
    energy = energy_corrections
  else
    allocate(energy(no_states,0), stat=ialloc); call err(ialloc)
  endif
  
  if (present(state_corrections)) then
    states = state_corrections
  else
    allocate(states(no_states,no_states,0), stat=ialloc); call err(ialloc)
  endif
  
  ! Check that input arrays are of consistent sizes.
  if ( size(perturbation,1) /= no_states .or. &
     & size(perturbation,2) /= no_states .or. &
     & size(energy,1)       /= no_states .or. &
     & size(states,1)       /= no_states .or. &
     & size(states,2)       /= no_states) then
    call err()
  endif
  
  if (size(energy,2)/=size(states,3)) then
    call err()
  endif
  
  ! Calculate the order of the next energy correction.
  order = size(energy,2)+1
  
  ! Allocate output.
  allocate(output(no_states), stat=ialloc); call err(ialloc)
  
  ! Calculate output.
  if (order==1) then
    do i=1,no_states
      output(i) = perturbation(i,i)
    enddo
  else
    do i=1,no_states
      output(i) = 0
      
      do j=1,no_states
        if (j/=i) then
          output(i) = output(i) + perturbation(i,j)*states(i,j,order-1)
        endif
      enddo
      
      do j=2,order-1
        output(i) = output(i) - energy(i,j)*states(i,i,order-j)
      enddo
    enddo
  endif
end function

! Calculates corrections to the states.
! Only returns the magnitude of the correction.
! The phase of the correction state(i,j) is the negative of the phase of
!    perturbation(i,j).
function calculate_state_correction(eigenvalues,perturbation, &
   & energy_corrections,state_corrections) result(output)
  implicit none
  
  real(dp), intent(in)           :: eigenvalues(:)
  real(dp), intent(in)           :: perturbation(:,:)
  real(dp), intent(in), optional :: energy_corrections(:,:)
  real(dp), intent(in), optional :: state_corrections(:,:,:)
  real(dp), allocatable          :: output(:,:)
  
  ! Previous corrections.
  real(dp), allocatable :: energy(:,:)
  real(dp), allocatable :: states(:,:,:)
  
  ! Array lengths.
  integer :: no_states
  integer :: order
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  no_states = size(eigenvalues)
  
  if (present(energy_corrections)) then
    energy = energy_corrections
  else
    allocate(energy(no_states,0), stat=ialloc); call err(ialloc)
  endif
  
  if (present(state_corrections)) then
    states = state_corrections
  else
    allocate(states(no_states,no_states,0), stat=ialloc); call err(ialloc)
  endif
  
  ! Check that input arrays are of consistent sizes.
  if ( size(perturbation,1) /= no_states .or. &
     & size(perturbation,2) /= no_states .or. &
     & size(energy,1)       /= no_states .or. &
     & size(states,1)       /= no_states .or. &
     & size(states,2)       /= no_states) then
    call err()
  endif
  
  if (size(energy,2)/=size(states,3)+1) then
    call err()
  endif
  
  ! Calculate the order of the next state correction.
  order = size(states,3)+1
  
  ! Allocate output.
  allocate(output(no_states,no_states), stat=ialloc); call err(ialloc)
  
  if (order==1) then
    do i=1,no_states
      output(i,i) = 0
      
      do j=1,no_states
        if (j/=i) then
          output(i,j) = perturbation(j,i)/(eigenvalues(i)-eigenvalues(j))
        endif
      enddo
    enddo
  else
    do i=1,no_states
      output(i,i) = 0
      do j=1,no_states
        do k=1,order-1
          output(i,i) = output(i,i) - 0.5_dp*states(i,j,k)*states(i,j,order-k)
        enddo
      enddo
      
      do j=1,no_states
        if (j/=i) then
          output(i,j) = 0
          
          do k=1,no_states
            if (k/=j) then
              output(i,j) = output(i,j) + perturbation(j,k)*states(i,k,order-1)
            endif
          enddo
          
          do k=2,order-1
            output(i,j) = output(i,j) - energy(i,k)*states(i,j,order-k)
          enddo
          
          output(i,j) = output(i,j)/(eigenvalues(i)-eigenvalues(j))
        endif
      enddo
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Calculates the correction due to a real perturbation.
! ----------------------------------------------------------------------
function calculate_perturbation_real(eigenvalues,perturbation,energy_order, &
   & state_order) result(output)
  implicit none
  
  real(dp), intent(in)   :: eigenvalues(:)
  real(dp), intent(in)   :: perturbation(:,:)
  integer,  intent(in)   :: energy_order
  integer,  intent(in)   :: state_order
  type(RealPerturbation) :: output
  
  integer               :: no_states
  real(dp), allocatable :: energy(:,:)
  real(dp), allocatable :: states(:,:,:)
  
  integer :: i,j,ialloc
  
  ! Check that the requested orders of perturbation are consistent.
  if (state_order/=energy_order .and. state_order/=energy_order-1) then
    call err()
  endif
  
  ! Check that the input array sizes are consistent.
  no_states = size(eigenvalues)
  if (any(shape(perturbation)/=no_states)) then
    call err()
  endif
  
  ! Allocate space for energy and state corrections.
  allocate( energy(no_states,energy_order),          &
          & states(no_states,no_states,state_order), &
          & stat=ialloc); call err(ialloc)
  
  ! Run perturbation theory.
  do i=1,energy_order
    energy(:,i) = calculate_energy_correction( eigenvalues,    &
                                             & perturbation,   &
                                             & energy(:,:i-1), &
                                             & states(:,:,:i-1))
    if (i<=state_order) then
      states(:,:,i) = calculate_state_correction( eigenvalues,  &
                                                & perturbation, &
                                                & energy(:,:i), &
                                                & states(:,:,:i-1))
    endif
  enddo
  
  do i=1,energy_order
    call print_line('')
    call print_line('Order '//i)
    call print_line('Energy corrections:')
    do j=1,no_states
      call print_line(energy(j,i))
    enddo
  enddo
  
  ! Generate output.
  allocate( output%energy(no_states),           &
          & output%states(no_states,no_states), &
          & stat=ialloc); call err(ialloc)
  output%energy = sum(energy,2)
  output%states = sum(states,3)
end function

! ----------------------------------------------------------------------
! Calculates the correction due to a complex perturbation.
! ----------------------------------------------------------------------
! The magnitude of the state correction is the same as the real case, but
!    the phase of each element is the negative of the equivalent element in
!    the perturbation.
function calculate_perturbation_complex(eigenvalues,perturbation, &
   & energy_order,state_order) result(output)
  implicit none
  
  real(dp),    intent(in)   :: eigenvalues(:)
  complex(dp), intent(in)   :: perturbation(:,:)
  integer,     intent(in)   :: energy_order
  integer,     intent(in)   :: state_order
  type(ComplexPerturbation) :: output
  
  integer                :: no_states
  type(RealPerturbation) :: real_output
  real(dp)               :: phase
  
  integer :: i,j,ialloc
  
  no_states = size(eigenvalues)
  
  real_output = calculate_perturbation( eigenvalues,       &
                                      & abs(perturbation), &
                                      & energy_order,      &
                                      & state_order)
  output%energy = real_output%energy
  allocate(output%states(no_states,no_states), stat=ialloc); call err(ialloc)
  do i=1,no_states
    do j=1,no_states
      phase = atan2(-aimag(perturbation(j,i)), real(perturbation(j,i)))
      output%states(j,i) = cmplx( real_output%states(j,i)*cos(phase), &
                                & real_output%states(j,i)*sin(phase), &
                                & dp)
    enddo
  enddo
end function
end module

! ======================================================================
! Calculates corrections due to perturbation theory.
! ======================================================================
module perturbation_module
  use constants_module, only : dp
  use string_module
  use io_module
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
end module

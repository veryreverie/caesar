! ======================================================================
! Calculates the frequencies of harmonic ground-states along each mode
!    which minimise the total energy of that ground-state.
! ======================================================================
module calculate_frequencies_module
  use common_module
  
  use states_module
  
  use potential_module
  use potential_pointer_module
  use anharmonic_data_module
  implicit none
  
  private
  
  public :: calculate_frequencies
  
  type :: FrequencyArray
    real(dp), allocatable :: frequencies(:)
  end type
contains

! The output is give as an array indexed as
!    anharmonic_data%degenerate_subspaces.
! (Every mode in a given subspace has the same frequency.)
subroutine calculate_frequencies(potential,anharmonic_data, &
   & frequency_convergence)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp),                 allocatable :: frequencies(:)
  
  type(SubspaceState), allocatable :: subspace_states(:)
  type(SubspaceState), allocatable :: states(:)
  
  type(PotentialPointer) :: new_potential
  
  type(FrequencyArray), allocatable :: old_frequencies(:)
  type(FrequencyArray), allocatable :: new_frequencies(:)
  
  real(dp) :: frequency
  
  integer :: i,j,k,ialloc
  
  subspaces = anharmonic_data%degenerate_subspaces
  
  ! Generate the first guess at the frequencies:
  !    - harmonic frequency             if > frequency_of_max_displacement.
  !    - frequency_of_max_displacement  otherwise.
  allocate( frequencies(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    frequencies(i) = max( subspaces(i)%frequency,                       &
                        & anharmonic_data%frequency_of_max_displacement )
  enddo
  
  ! Generate ground states at first guess frequencies.
  allocate(states(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    subspace_states = generate_subspace_states( &
               & subspaces(i),                  &
               & frequencies(i),                &
               & anharmonic_data%complex_modes, &
               & 0                              )
    if (size(subspace_states)/=1) then
      call err()
    endif
    states(i) = subspace_states(1)
  enddo
  
  old_frequencies = [FrequencyArray(frequencies)]
  new_frequencies = [FrequencyArray::]
  i = 1
  do
    ! Append an empty array to new_frequencies.
    new_frequencies = [ new_frequencies,                               &
                      & FrequencyArray([(0.0_dp,j=1,size(subspaces))]) ]
    
    ! Calculate update frequencies.
    do j=1,size(subspaces)
      ! Integrate the potential across all other subspaces.
      new_potential = potential
      do k=1,size(subspaces)
        if (k/=j) then
          call new_potential%braket(states(k),states(k),anharmonic_data)
        endif
      enddo
      
      ! Set the constant energy to zero, to stabilise minima finding.
      call new_potential%zero_energy()
      
      ! Find the frequency which minimises total energy.
      new_frequencies(i)%frequencies(j) = optimise_frequency( &
                                      & new_potential,        &
                                      & states(j),            &
                                      & anharmonic_data,      &
                                      & frequency_convergence )
    enddo
    
    ! Use a Pulay scheme to update frequencies self-consistently.
    old_frequencies = [ old_frequencies,                       &
                      & pulay(old_frequencies,new_frequencies) ]
    
    ! Check whether the Pulay iteration has converged.
    i = i+1
    if (sum(abs( old_frequencies(i)%frequencies &
             & - old_frequencies(i-1)%frequencies))<frequency_convergence) then
      frequencies = old_frequencies(i)%frequencies
      exit
    endif
  enddo
  
  ! TODO
  ! Generate basis states from frequencies.
  ! Return basis.
  
end subroutine

! Find the frequency which minimises energy in a single subspace.
recursive function optimise_frequency(potential,state,anharmonic_data, &
   & frequency_convergence) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(SubspaceState),  intent(in) :: state
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  real(dp)                         :: output
  
  real(dp) :: frequencies(3)
  real(dp) :: energies(3)
  
  type(SubspaceState)    :: new_state
  type(PotentialPointer) :: new_potential
  
  real(dp) :: first_derivative
  real(dp) :: second_derivative
  
  real(dp) :: old_frequency
  real(dp) :: new_frequency
  
  integer :: i
  
  old_frequency   = state%frequency
  new_state       = state
  
  ! Calculate [w-dw, w, w+dw].
  frequencies = old_frequency
  frequencies(1) = frequencies(1) - 0.01_dp*frequency_convergence
  frequencies(3) = frequencies(3) + 0.01_dp*frequency_convergence
  
  ! Calculate [U(w-dw), U(w), U(w+dw)].
  do i=1,3
    new_potential = potential
    new_state%frequency = frequencies(i)
    call new_potential%braket(new_state,new_state,anharmonic_data)
    energies(i) = new_potential%undisplaced_energy()
  enddo
  
  ! Calculate dU/dw = (U(w+dw)-U(w-dw))/(2dw).
  first_derivative = (energies(3)-energies(1)) &
                 & / (0.02_dp*frequency_convergence)
  
  ! Calculate d2U/dw2 = (U(w+dw)+U(w-dw)-2U(w))/(dw)^2.
  second_derivative = (energies(1)+energies(3)-2*energies(2)) &
                  & / (0.01_dp*frequency_convergence)**2
  
  ! At the extrema, (w=w1), dU/dw=0, so w1 = w - (dU/dw)/(d2U/dw2).
  ! If |w1-w|>w/2, or if dU/dw<0 then cap |w1-w| at w/2.
  if (abs(old_frequency)*second_derivative<=abs(first_derivative)) then
    if (first_derivative>0) then
      new_frequency = 0.5_dp * old_frequency
    elseif (first_derivative<0) then
      new_frequency = 1.5_dp * old_frequency
    else
      new_frequency = old_frequency
      output = new_frequency
      return
    endif
    
    new_state%frequency = new_frequency
    output = optimise_frequency( potential,            &
                               & new_state,            &
                               & anharmonic_data,      &
                               & frequency_convergence )
  else
    new_frequency = old_frequency - first_derivative/second_derivative
    if (abs(new_frequency-old_frequency)<frequency_convergence) then
      output = new_frequency
      return
    else
      new_state%frequency = new_frequency
      output = optimise_frequency( potential,            &
                                 & new_state,            &
                                 & anharmonic_data,      &
                                 & frequency_convergence )
    endif
  endif
end function

! Use a Pulay scheme to minimise frequencies of all subspaces.
function pulay(old_frequencies,new_frequencies) result(output)
  implicit none
  
  type(FrequencyArray), intent(in) :: old_frequencies(:)
  type(FrequencyArray), intent(in) :: new_frequencies(:)
  type(FrequencyArray)             :: output
  
  output = old_frequencies(size(old_frequencies))
  
  ! TODO: Write Pulay scheme.
end function
end module

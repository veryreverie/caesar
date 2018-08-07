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
  
  type(RealVector), allocatable :: old_frequencies(:)
  type(RealVector), allocatable :: new_frequencies(:)
  
  real(dp) :: frequency
  
  integer :: first_pulay_iteration
  
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
  
  ! Find self-consistent frequencies which minimise the energy.
  old_frequencies = [vec(frequencies)]
  new_frequencies = [RealVector::]
  i = 1
  iter: do
    ! Update the states to have the new frequencies.
    call print_line('')
    call print_line('========================================')
    call print_line(old_frequencies(i))
    states%frequency = dble(old_frequencies(i))
    call print_line(states%frequency)
    
    ! Calculate mean-field potentials from the old frequencies, and use these
    !    mean-field potentials to calculate new frequencies.
    new_frequencies = [ new_frequencies,                              &
                      & optimise_frequencies( potential,              &
                      &                       states,                 &
                      &                       anharmonic_data,        &
                      &                       frequency_convergence ) ]
    call print_line(new_frequencies(i))
    
    ! Use a Pulay scheme to converge towards the self-consistent solution,
    !    such that new frequencies = old frequencies.
    if (i<=200) then
      old_frequencies = [ old_frequencies,   &
                        & new_frequencies(i) ]
    else
      old_frequencies = [ old_frequencies,                                &
                        & pulay(old_frequencies(2:), new_frequencies(2:)) ]
    endif
    
    call print_line(old_frequencies(i+1))
    call print_line('========================================')
    
    ! Check whether the frequencies have converged.
    if (i>2) then
      if ( sum(abs(dble(old_frequencies(i)-old_frequencies(i-1)))) &
       & < frequency_convergence                                   ) then
        frequencies = dble(old_frequencies(i))
        exit iter
      endif
    endif
    
    ! Increment the loop counter.
    i = i+1
  enddo iter
  
  ! TODO
  ! Generate basis states from frequencies.
  ! Return basis.
  
end subroutine

! For each subspace, integrate the potential across all other subspaces
!    to get a single-subspace mean-field potential.
! Then find the subspace frequency which minimises the single-subspace energy.
recursive function optimise_frequencies(potential,states,anharmonic_data, &
   & frequency_convergence) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(SubspaceState),  intent(in) :: states(:)
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  type(RealVector)                 :: output
  
  real(dp), allocatable  :: new_frequencies(:)
  type(PotentialPointer) :: new_potential
  
  integer :: i,j,ialloc
  
  new_frequencies = [(0.0_dp,i=1,size(states))]
  
  ! Calculate update frequencies.
  do i=1,size(states)
    ! Integrate the potential across all other subspaces.
    new_potential = potential
    do j=1,size(states)
      if (j/=i) then
        call new_potential%braket(states(j),states(j),anharmonic_data)
      endif
    enddo
    
    ! Set the constant energy to zero, to stabilise minima finding.
    call new_potential%zero_energy()
    call print_lines(new_potential)
    
    ! Find the frequency which minimises total energy.
    new_frequencies(i) = optimise_frequency( new_potential,        &
                                           & states(i),            &
                                           & anharmonic_data,      &
                                           & frequency_convergence )
  enddo
  
  output = new_frequencies
end function

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
  
  old_frequency = state%frequency
  new_state     = state
  
  ! Calculate [w-dw, w, w+dw].
  frequencies = old_frequency
  frequencies(1) = frequencies(1) - 0.01_dp*frequency_convergence
  frequencies(3) = frequencies(3) + 0.01_dp*frequency_convergence
  
  ! Calculate [U(w-dw), U(w), U(w+dw)].
  do i=1,3
    new_potential = potential
    new_state%frequency = frequencies(i)
    call new_potential%braket(new_state,new_state,anharmonic_data)
    energies(i) = new_potential%undisplaced_energy() &
              & + kinetic_energy(new_state, new_state)
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
end module

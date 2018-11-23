! ======================================================================
! Generates a basis of states which span each degenerate subspace.
! ======================================================================
! The states are products single- and double-mode states.
! The single-mode states are in two forms:
!    - monomial states: (u_i)^(n_i)|0>.
!    - harmonic states: (a'_i)^(n_i)|0>.
! The double-mode states are in two forms:
!    - monomial states: (u_i)^(n_i)(u_j)^(n_j)|0,0>.
!    - harmonic states: (a'_i)^(n_i)(a'_j)^(n_j)|0,0>.
! In both cases, the states |0> and |0,0> are the ground states of an
!    effective harmonic potential with harmonic frequencies chosen such that
!    the energy of |0> and |0,0>  w.r.t. the anharmonic potential is minimised.
! N.B. the monomial states are normalised, but they are NOT orthogonal.
module generate_basis_module
  use common_module
  
  use states_module
  
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: generate_basis
contains

function generate_basis(potential,anharmonic_data,frequency_convergence, &
   & max_pulay_iterations,pre_pulay_iterations,pre_pulay_damping,        &
   & no_basis_states) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  integer,              intent(in) :: max_pulay_iterations
  integer,              intent(in) :: pre_pulay_iterations
  real(dp),             intent(in) :: pre_pulay_damping
  integer,              intent(in) :: no_basis_states
  type(SubspaceBasis), allocatable :: output(:)
  
  type(DegenerateSubspace), allocatable :: subspaces(:)
  real(dp),                 allocatable :: frequencies(:)
  
  type(MonomialState), allocatable :: states(:)
  type(MonomialState), allocatable :: subspace_states(:)
  
  type(RealVector), allocatable :: old_frequencies(:)
  type(RealVector), allocatable :: new_frequencies(:)
  
  type(RealVector) :: next_old_frequencies
  
  integer :: i,j,ialloc
  
  if (pre_pulay_damping<0_dp .or. pre_pulay_damping>1_dp) then
    call print_line(ERROR//': pre_pulay_damping must be between 0 and 1.')
    stop
  elseif (pre_pulay_iterations<0) then
    call print_line(ERROR//': pre_pulay_iterations must be >= 0.')
    stop
  elseif (max_pulay_iterations<0) then
    call print_line(ERROR//': max_pulay_iterations must be >= 0.')
    stop
  endif
  
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
    subspace_states = generate_monomial_states( &
               & subspaces(i),                  &
               & frequencies(i),                &
               & anharmonic_data%complex_modes, &
               & maximum_power = 0              )
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
    states%frequency = dble(old_frequencies(i))
    
    ! Calculate mean-field potentials from the old frequencies, and use these
    !    mean-field potentials to calculate new frequencies.
    new_frequencies = [ new_frequencies,                              &
                      & optimise_frequencies( potential,              &
                      &                       states,                 &
                      &                       anharmonic_data,        &
                      &                       frequency_convergence ) ]
    
    ! Use a damped iterative scheme or a Pulay scheme to converge towards the
    !    self-consistent solution where new frequencies = old frequencies.
    if (i<=pre_pulay_iterations) then
      next_old_frequencies = (1-pre_pulay_damping) * old_frequencies(i) &
                         & + pre_pulay_damping * new_frequencies(i)
    else
      j = max(2, i-max_pulay_iterations+1)
      next_old_frequencies = pulay(old_frequencies(j:), new_frequencies(j:))
    endif
    
    old_frequencies = [old_frequencies, next_old_frequencies]
    
    ! Increment the loop counter.
    i = i+1
    
    ! Check whether the frequencies have converged.
    ! If the frequencies aren't converged, return to the top of the loop.
    if (i<=5) then
      cycle iter
    else
      do j=1,5
        if ( sum(abs(dble(old_frequencies(i)-old_frequencies(i-j)))) &
         & > frequency_convergence                                   ) then
          cycle iter
        endif
      enddo
    endif
    
    ! If the frequencies are converged, break out of the loop.
    frequencies = dble(old_frequencies(i))
    exit iter
  enddo iter
  
  ! Use calculated effective harmonic frequencies to generate bases.
  allocate(output(size(subspaces)), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    output(i) = SubspaceBasis( subspaces(i),                              &
                             & frequencies(i),                            &
                             & anharmonic_data%complex_modes,             &
                             & anharmonic_data%qpoints,                   &
                             & anharmonic_data%anharmonic_supercell,      &
                             & no_basis_states-1,                         &
                             & anharmonic_data%potential_expansion_order, &
                             & anharmonic_data%structure%symmetries       )
  enddo
end function

! For each subspace, integrate the potential across all other subspaces
!    to get a single-subspace mean-field potential.
! Then find the subspace frequency which minimises the single-subspace energy.
function optimise_frequencies(potential,states,anharmonic_data, &
   & frequency_convergence) result(output)
  implicit none
  
  class(PotentialData), intent(in) :: potential
  type(MonomialState),  intent(in) :: states(:)
  type(AnharmonicData), intent(in) :: anharmonic_data
  real(dp),             intent(in) :: frequency_convergence
  type(RealVector)                 :: output
  
  type(PotentialPointer), allocatable :: subspace_potentials(:)
  
  real(dp), allocatable  :: new_frequencies(:)
  
  integer :: i
  
  ! Calculate the single-subspace potentials, defined as
  !    V_i = (prod_{j/=i}<j|)V(prod_{j/=i}|j>).
  subspace_potentials = generate_subspace_potentials( potential,      &
                                                    & states,         &
                                                    & anharmonic_data )
  
  new_frequencies = [(0.0_dp,i=1,size(states))]
  
  ! Calculate update frequencies.
  do i=1,size(states)
    ! Find the frequency which minimises total energy.
    new_frequencies(i) = optimise_frequency(      &
       & subspace_potentials(i),                  &
       & states(i),                               &
       & anharmonic_data,                         &
       & anharmonic_data%degenerate_subspaces(i), &
       & frequency_convergence )
  enddo
  
  output = new_frequencies
end function

! Find the frequency which minimises energy in a single subspace.
function optimise_frequency(potential,state,anharmonic_data,subspace, &
   & frequency_convergence) result(output)
  implicit none
  
  class(PotentialData),     intent(in) :: potential
  type(MonomialState),      intent(in) :: state
  type(AnharmonicData),     intent(in) :: anharmonic_data
  type(DegenerateSubspace), intent(in) :: subspace
  real(dp),                 intent(in) :: frequency_convergence
  real(dp)                             :: output
  
  real(dp) :: frequencies(3)
  real(dp) :: energies(3)
  
  type(MonomialState)    :: new_state
  
  real(dp) :: first_derivative
  real(dp) :: second_derivative
  
  real(dp) :: old_frequency
  real(dp) :: new_frequency
  
  integer :: i
  
  old_frequency = state%frequency
  new_state = state
  
  do
    ! Calculate [w-dw, w, w+dw].
    frequencies = [ old_frequency - 0.01_dp*frequency_convergence, &
                  & old_frequency,                                 &
                  & old_frequency + 0.01_dp*frequency_convergence  ]
    
    ! Calculate [U(w-dw), U(w), U(w+dw)].
    do i=1,3
      new_state%frequency = frequencies(i)
      energies(i) = potential%potential_energy( new_state,               &
                &                               new_state,               &
                &                               anharmonic_data )        &
                & + kinetic_energy( new_state,                           &
                &                   new_state,                           &
                &                   subspace,                            &
                &                   anharmonic_data%anharmonic_supercell )
    enddo
    
    ! Calculate dU/dw = (U(w+dw)-U(w-dw))/(2dw).
    first_derivative = (energies(3)-energies(1)) &
                   & / (0.02_dp*frequency_convergence)
    
    ! Calculate d2U/dw2 = (U(w+dw)+U(w-dw)-2U(w))/(dw)^2.
    second_derivative = (energies(1)+energies(3)-2*energies(2)) &
                    & / (0.01_dp*frequency_convergence)**2
    
    ! Update the frequency, and check for convergence.
    ! At the extrema (w=w1), dU/dw=0. As such, w1 = w - (dU/dw)/(d2U/dw2).
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
    else
      new_frequency = old_frequency - first_derivative/second_derivative
      if (abs(new_frequency-old_frequency)<frequency_convergence) then
        output = new_frequency
        return
      else
      endif
    endif
    
    ! If the frequency hasn't converged, update the frequency,
    !    and return to the start of the loop.
    old_frequency = new_frequency
  enddo
end function
end module

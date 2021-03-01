submodule (caesar_generate_subspace_potentials_module) caesar_generate_subspace_potentials_submodule
  use caesar_anharmonic_module
contains

module procedure generate_subspace_potentials
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
  ! The correction energy, -<V>(n-1)/n
  type(PotentialPointer) :: integrated_potential
  real(dp)               :: correction_energy
  
  integer :: i,j,ialloc
  
  if (size(subspace_states)==0) then
    call print_line(ERROR//': No states to integrate over.')
    call err()
  endif
  
  ! Initialise to a single interval spanning all states,
  !    and the un-integrated potential.
  mins_in = [1]
  maxs_in = [size(subspace_states)]
  allocate(output(size(subspace_states)), stat=ialloc); call err(ialloc)
  output(1) = PotentialPointer(potential)
  call output(1)%zero_energy()
  
  ! Loop over iteration until every interval contains exactly one subspace.
  do while (any(mins_in/=maxs_in))
    mins_out = [integer::]
    maxs_out = [integer::]
    
    ! Loop over intervals.
    do i=1,size(mins_in)
      if (mins_in(i)==maxs_in(i)) then
        ! The interval contains only one subspace; nothing needs doing.
        mins_out = [mins_out, mins_in(i)]
        maxs_out = [maxs_out, maxs_in(i)]
      else
        ! The interval contains more than one subspace;
        !    and generate two integrated potentials.
        
        ! Bisect the interval (or split it as evenly as possible).
        s = mins_in(i) + (maxs_in(i)-mins_in(i)+1)/2
        mins_out = [mins_out, mins_in(i), s     ]
        maxs_out = [maxs_out, s-1   , maxs_in(i)]
        
        ! Copy the potential from the first interval to the second.
        output(s) = output(mins_in(i))
        
        ! Integrate the first potential over all states in the second interval,
        !    and the second potential over all states in the first interval.
        do j=s,maxs_in(i)
          call output(mins_in(i))%braket(            &
             & states          = subspace_states(j), &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
        do j=mins_in(i),s-1
          call output(s)%braket(                     &
             & states          = subspace_states(j), &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Calculate the fully-integrated potential,
  !    <V> = (prod_i<i|)V(prod_i|i>).
  ! Use this to correct for the overcounting of terms.
  integrated_potential = output(1)
  
  call integrated_potential%braket( states          = subspace_states(1), &
                                  & subspace        = subspaces(1),       &
                                  & subspace_basis  = subspace_bases(1),  &
                                  & anharmonic_data = anharmonic_data     )
  
  correction_energy = integrated_potential%undisplaced_energy() &
                  & * (1.0_dp-size(subspaces))/size(subspaces)
  
  call output%add_constant(correction_energy)
  
  ! Process the subspace potentials if necessary.
  call subspace_bases%process_subspace_potential( output,          &
                                                & subspace_states, &
                                                & subspaces,       &
                                                & anharmonic_data  )
  
  do i=1,size(output)
    if (present(old_subspace_potentials)) then
      call output(i)%optimise_subspace_potential( subspaces(i),               &
                                                & subspace_bases(i),          &
                                                & old_subspace_potentials(i), &
                                                & anharmonic_data             )
    else
      call output(i)%optimise_subspace_potential( &
              & subspaces(i),                     &
              & subspace_bases(i),                &
              & anharmonic_data = anharmonic_data )
    endif
  enddo
end procedure

module procedure generate_subspace_stresses
  ! The minimum and maximum indices in each interval.
  integer, allocatable :: mins_in(:)
  integer, allocatable :: maxs_in(:)
  integer, allocatable :: mins_out(:)
  integer, allocatable :: maxs_out(:)
  
  ! The half-way point of each interval.
  integer :: s
  
  ! The correction stress, -<S>(n-1)/n
  type(StressPointer) :: integrated_stress
  type(RealMatrix)    :: correction_stress
  
  integer :: i,j,ialloc
  
  if (size(subspace_states)==0) then
    call print_line(ERROR//': No states to integrate over.')
    call err()
  endif
  
  ! Initialise to a single interval spanning all states,
  !    and the un-integrated stress.
  mins_in = [1]
  maxs_in = [size(subspace_states)]
  allocate(output(size(subspace_states)), stat=ialloc); call err(ialloc)
  output(1) = StressPointer(stress)
  call output(1)%zero_stress()
  
  ! Loop over iteration until every interval contains exactly one subspace.
  do while (any(mins_in/=maxs_in))
    mins_out = [integer::]
    maxs_out = [integer::]
    
    ! Loop over intervals.
    do i=1,size(mins_in)
      if (mins_in(i)==maxs_in(i)) then
        ! The interval contains only one subspace; nothing needs doing.
        mins_out = [mins_out, mins_in(i)]
        maxs_out = [maxs_out, maxs_in(i)]
      else
        ! The interval contains more than one subspace;
        !    and generate two integrated stresses.
        
        ! Bisect the interval (or split it as evenly as possible).
        s = mins_in(i) + (maxs_in(i)-mins_in(i)+1)/2
        mins_out = [mins_out, mins_in(i), s     ]
        maxs_out = [maxs_out, s-1   , maxs_in(i)]
        
        ! Copy the stress from the first interval to the second.
        output(s) = output(mins_in(i))
        
        ! Integrate the first stress over all states in the second interval,
        !    and the second stress over all states in the first interval.
        do j=s,maxs_in(i)
          call output(mins_in(i))%braket(            &
             & states          = subspace_states(j), &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
        do j=mins_in(i),s-1
          call output(s)%braket(                     &
             & states          = subspace_states(j), &
             & subspace        = subspaces(j),       &
             & subspace_basis  = subspace_bases(j),  &
             & anharmonic_data = anharmonic_data     )
        enddo
      endif
    enddo
    
    mins_in = mins_out
    maxs_in = maxs_out
  enddo
  
  ! Calculate the fully-integrated potential,
  !    <V> = (prod_i<i|)V(prod_i|i>).
  ! Use this to correct for the overcounting of terms.
  integrated_stress = output(1)
  
  call integrated_stress%braket( states          = subspace_states(1), &
                               & subspace        = subspaces(1),       &
                               & subspace_basis  = subspace_bases(1),  &
                               & anharmonic_data = anharmonic_data     )
  
  correction_stress = integrated_stress%undisplaced_stress() &
                  & * (1.0_dp-size(subspaces))/size(subspaces)
  
  call output%add_constant(correction_stress)
  
  ! Process the subspace stresses if necessary.
  call subspace_bases%process_subspace_stress( output,          &
                                             & subspace_states, &
                                             & subspaces,       &
                                             & anharmonic_data  )
end procedure
end submodule

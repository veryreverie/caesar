submodule (caesar_vscf_module) caesar_vscf_submodule
  use caesar_anharmonic_module
contains

module procedure new_VscfOutput
  this%potential = PotentialPointer(potential)
  this%states = BasisStatesPointer(states)
  if (present(stress)) then
    this%stress = StressPointer(stress)
  endif
end procedure

module procedure run_vscf
  ! Single-subspace potentials and states.
  type(PotentialPointer),   allocatable :: subspace_potentials(:)
  type(BasisStatesPointer), allocatable :: subspace_states(:)
  type(StressPointer),      allocatable :: subspace_stresses(:)
  
  ! Potential coefficients.
  integer,                    allocatable :: first_coefficient(:)
  integer,                    allocatable :: last_coefficient(:)
  real(dp),                   allocatable :: coefficients(:)
  real(dp)                                :: free_energy
  type(PotentialBasePointer), allocatable :: variable_basis_functions(:)
  real(dp),                   allocatable :: free_energy_gradient(:)
  
  ! The Pulay scheme solver.
  type(PulaySolver) :: solver
  
  ! Thermodynamic variables.
  type(ThermodynamicData), allocatable :: thermodynamic_data(:)
  real(dp),                allocatable :: free_energies(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! Generate initial single-subspace states,
  !    and use these states to generate initial single-subpsace potentials.
  if (.not. present(starting_configuration)) then
    call print_line('Generating initial states and potentials.')
    subspace_states = BasisStatesPointer(subspace_bases%initial_states( &
                                                      & subspaces,      &
                                                      & thermal_energy, &
                                                      & anharmonic_data ))
    subspace_potentials = generate_subspace_potentials( &
                    & potential,                        &
                    & subspaces,                        &
                    & subspace_bases,                   &
                    & subspace_states,                  &
                    & anharmonic_data = anharmonic_data )
  else
    call print_line('Using previous configuration for initial states and &
       &potentials.')
    subspace_states = starting_configuration%states
    subspace_potentials = starting_configuration%potential
  endif
  
  ! Construct a single vector of coefficients from the subspace potentials.
  allocate( coefficients(0),                    &
          & first_coefficient(size(subspaces)), &
          & last_coefficient(size(subspaces)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    first_coefficient(i) = size(coefficients)+1
    coefficients = [ coefficients,                                        &
                   & subspace_potentials(i)%coefficients(anharmonic_data) ]
    last_coefficient(i) = size(coefficients)
  enddo
  allocate( free_energy_gradient(size(coefficients)), &
          & stat=ialloc); call err(ialloc)
  
  ! Write initial potentials to file.
  if (present(convergence_file)) then
    call convergence_file%print_line('Convergence Threshold: '//          &
                                    & convergence_data%energy_convergence )
    call convergence_file%print_line(repeat('=',50))
    call convergence_file%print_line('Initial Subspace Potentials:')
    call convergence_file%print_line(repeat('=',50))
    do i=1,size(subspace_potentials)
      call convergence_file%print_line(repeat('-',50))
      call convergence_file%print_line('Subspace '//subspaces(i)%id)
      call convergence_file%print_line('Subspace coefficients: '// &
         & first_coefficient(i)//' to '//last_coefficient(i))
      call convergence_file%print_line(repeat('-',50))
      call convergence_file%print_lines(subspace_potentials(i))
      call convergence_file%print_line(repeat('-',50))
    enddo
    call convergence_file%print_line(repeat('=',50))
  endif
  
  ! Initialise Pulay solver.
  solver = PulaySolver( convergence_data, &
                      & random_generator, &
                      & coefficients      )
  
  ! Run Pulay scheme.
  i = 1
  free_energies = [real(dp)::]
  do
    ! Generate single-subspace potentials from Pulay scheme.
    call solver%calculate_x()
    coefficients = solver%get_x()
    do j=1,size(subspaces)
      call subspace_potentials(j)%set_coefficients(                &
         & coefficients(first_coefficient(j):last_coefficient(j)), &
         & anharmonic_data                                         )
    enddo
    
    ! Write coefficients to file.
    if (present(convergence_file)) then
      call convergence_file%print_line('Iteration '//i)
      call convergence_file%print_line('Input coefficients:')
      call convergence_file%print_line(coefficients)
    endif
    
    ! Use single-subspace potentials to calculate new single-subspace states.
    call print_line('Generating single-subspace ground states.')
    subspace_states =                                              &
       & subspace_bases%calculate_states( subspaces,               &
       &                                  subspace_potentials,     &
       &                                  thermal_energy,          &
       &                                  state_energy_cutoff,     &
       &                                  convergence_data,        &
       &                                  anharmonic_data          )
    
    ! Use single-subspace states to calculate new single-subspace potentials.
    call print_line('Generating single-subspace potentials.')
    subspace_potentials = generate_subspace_potentials( potential,           &
                                                      & subspaces,           &
                                                      & subspace_bases,      &
                                                      & subspace_states,     &
                                                      & subspace_potentials, &
                                                      & anharmonic_data      )
    
    do j=1,size(subspaces)
      coefficients(first_coefficient(j):last_coefficient(j)) = &
         & subspace_potentials(j)%coefficients(anharmonic_data)
    enddo
    
    ! Calculate the free energy from the potentials and states.
    thermodynamic_data = subspace_bases%thermodynamic_data( &
                          & thermal_energy,                 &
                          & subspace_states,                &
                          & subspaces,                      &
                          & subspace_potentials,            &
                          & anharmonic_data=anharmonic_data )
    free_energy = sum(thermodynamic_data%free_energy) &
              & / anharmonic_data%anharmonic_supercell%sc_size
    
    ! Update the Pulay scheme.
    call solver%set_f( coefficients,         &
                     & free_energy,          &
                     & print_progress=.true. )
    if (solver%gradient_requested()) then
      call print_line('Calculating free energy derivative.')
      do j=1,size(subspaces)
        variable_basis_functions = &
           & subspace_potentials(j)%variable_basis_functions(anharmonic_data)
        free_energy_gradient(first_coefficient(j):last_coefficient(j)) = &
           & subspace_bases(j)%free_energy_gradient(                     &
           &               subspace_potentials(j),                       &
           &               variable_basis_functions,                     &
           &               subspaces(j),                                 &
           &               subspace_states(j),                           &
           &               thermal_energy,                               &
           &               state_energy_cutoff,                          &
           &               anharmonic_data           )
      enddo
      call solver%set_gradient(free_energy_gradient)
    endif
    
    ! Write coefficients to file.
    if (present(convergence_file)) then
      call convergence_file%print_line('Output coefficients:')
      call convergence_file%print_line(coefficients)
      call convergence_file%print_line('Free energy:')
      call convergence_file%print_line(free_energy)
    endif
    
    ! Check for convergence.
    if (solver%converged()) then
      call print_line('Convergence reached.')
      
      ! Calculate subspace stresses.
      if (present(stress)) then
        subspace_stresses = generate_subspace_stresses( stress,          &
                                                      & subspaces,       &
                                                      & subspace_bases,  &
                                                      & subspace_states, &
                                                      & anharmonic_data  )
      endif
      
      ! Construct output and return.
      if (present(stress)) then
        output = VscfOutput( subspace_potentials, &
                           & subspace_states,     &
                           & subspace_stresses    )
      else
        output = VscfOutput( subspace_potentials, &
                           & subspace_states      )
      endif
      return
    endif
    
    ! Increment the loop counter.
    i = i+1
  enddo
end procedure
end submodule

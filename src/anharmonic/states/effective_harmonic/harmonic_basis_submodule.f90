submodule (caesar_harmonic_basis_module) caesar_harmonic_basis_submodule
  use caesar_effective_harmonic_module
contains

module procedure startup_harmonic_basis
  type(HarmonicBasis) :: basis
  
  call basis%startup()
end procedure

module procedure representation_HarmonicBasis
  output = 'harmonic'
end procedure

module procedure new_HarmonicBasis
  this%subspace_id = subspace_id
  this%frequency = frequency
  this%supercell_size = supercell_size
end procedure

module procedure initial_states_HarmonicBasis
  output = BasisStatesPointer(HarmonicStates( subspace%id,    &
                                            & this%frequency, &
                                            & thermal_energy  ))
end procedure

module procedure calculate_states_HarmonicBasis
  real(dp) :: energy_convergence
  
  type(NewtonRaphson) :: solver
  
  real(dp)                :: frequencies(3)
  type(ThermodynamicData) :: observables(3)
  
  real(dp) :: frequency
  
  integer :: i
  
  energy_convergence = convergence_data%energy_convergence
  
  solver = NewtonRaphson(                                  &
     !& starting_value        = this%frequency,             &
     & starting_value        = 1e-2_dp,                    &
     & finite_displacement   = 0.01_dp*energy_convergence, &
     & convergence_threshold = 0.5_dp*energy_convergence,  &
     & lower_bound           = 1e-300_dp                   )
  i = 0
  do 
    frequencies = solver%get_inputs()
    
    if (frequencies(2)<2e-300_dp) then
      call print_line(ERROR//': VSCHA frequency in subspace '//subspace%id// &
         & ' underflowed. This likely means the subspace potential is not &
         &well bounded from below.')
      call print_line('Subspace potential:')
      call print_lines(subspace_potential)
      call print_line(frequencies)
      call err()
    elseif (.not. all(abs(frequencies)<1e300_dp)) then
      call print_line(ERROR//': Newton-Raphson scheme diverged.')
      call print_line('Iteration   : '//i)
      call print_line('Frequency   : '//frequencies)
      call print_line('Free energy : '//observables%free_energy)
      call print_line('')
      call print_line('Subspace potential:')
      call print_lines(subspace_potential)
      call err()
    endif
    
    observables = effective_harmonic_observables( &
         & thermal_energy  = thermal_energy,      &
         & potential       = subspace_potential,  &
         & frequency       = frequencies,         &
         & num_dimensions  = size(subspace),      &
         & supercell_size  = this%supercell_size, &
         & anharmonic_data = anharmonic_data      )
    
    call solver%set_outputs(observables%free_energy)
    
    if (solver%converged()) then
      frequency = solver%solution()
      exit
    endif
    
    i = i+1
    if (modulo(i,1000)==0) then
      call print_line(WARNING//': Newton-Raphson scheme taking a long time to &
         &converge.')
      call print_line('Iteration   : '//i)
      call print_line('Frequency   : '//frequencies)
      call print_line('Free energy : '//observables%free_energy)
      call quit()
    endif
  enddo
  
  frequency = max(frequency, 1e-3_dp)
  
  output = BasisStatesPointer(HarmonicStates( this%subspace_id, &
                                            & frequency,        &
                                            & thermal_energy    ))
end procedure

module procedure mode_ids_HarmonicBasis
  type(ComplexMode), allocatable :: modes(:)
  
  modes = subspace%modes(anharmonic_data%complex_modes)
  
  modes = modes(filter(modes%id<=modes%paired_id))
  
  output = modes%id
end procedure

module procedure paired_mode_ids_HarmonicBasis
  type(ComplexMode), allocatable :: modes(:)
  
  modes = subspace%modes(anharmonic_data%complex_modes)
  
  modes = modes(filter(modes%id<=modes%paired_id))
  
  output = modes%paired_id
end procedure

module procedure inner_product_HarmonicBasis
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end procedure

module procedure integrate_BasisState_HarmonicBasis
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end procedure

module procedure kinetic_energy_HarmonicBasis
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end procedure

module procedure harmonic_potential_energy_HarmonicBasis
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end procedure

module procedure kinetic_stress_HarmonicBasis
  call print_line(CODE_ERROR//': Procedures involving individual states have &
     &not been implemented for HarmonicBasis.')
  call err()
end procedure

module procedure thermodynamic_data_HarmonicBasis
  type(HarmonicStates) :: harmonic_states
  
  harmonic_states = HarmonicStates(states)
  
  output = effective_harmonic_observables(            &
     & thermal_energy  = thermal_energy,              &
     & potential       = subspace_potential,          &
     & frequency       = harmonic_states%frequency,   &
     & num_dimensions  = size(subspace),              &
     & supercell_size  = this%supercell_size,         &
     & anharmonic_data = anharmonic_data              )
end procedure

module procedure wavefunctions_HarmonicBasis
  call err()
end procedure

module procedure integrate_BasisStates_HarmonicBasis
  type(HarmonicStates) :: harmonic_states
  
  harmonic_states = HarmonicStates(states)
  
  output = product(monomial%modes%harmonic_expectation( &
                      & harmonic_states%frequency,      &
                      & harmonic_states%thermal_energy, &
                      & this%supercell_size             ))
end procedure

module procedure free_energy_gradient_HarmonicBasis
  integer :: i
  
  ! TODO: implement this.
  output = [(0.0_dp,i=1,size(basis_functions))]
end procedure

module procedure read_HarmonicBasis
  integer  :: subspace_id
  real(dp) :: frequency
  integer  :: supercell_size
  
  select type(this); type is(HarmonicBasis)
    subspace_id = int(token(input(1),3))
    frequency = dble(token(input(2),3))
    supercell_size = int(token(input(3),4))
    
    this = HarmonicBasis(subspace_id, frequency, supercell_size)
  class default
    call err()
  end select
end procedure

module procedure write_HarmonicBasis
  select type(this); type is(HarmonicBasis)
    output = [ 'Subspace       : '//this%subspace_id,   &
             & 'Frequency      : '//this%frequency,     &
             & 'Supercell Size : '//this%supercell_size ]
  class default
    call err()
  end select
end procedure

module procedure new_HarmonicBasis_Strings
  call this%read(input)
end procedure

module procedure new_HarmonicBasis_StringArray
  this = HarmonicBasis(str(input))
end procedure
end submodule

submodule (caesar_full_subspace_basis_module) caesar_full_subspace_basis_submodule
  use caesar_states_module
contains

module procedure new_FullSubspaceBasis
  this%supercell_size  = supercell_size
  this%maximum_power   = maximum_power
  this%expansion_order = expansion_order
  this%subspace_id     = subspace_id
  this%frequency       = frequency
  this%wavevectors     = wavevectors
end procedure

module procedure representation_FullSubspaceBasis
  output = 'full_subspace'
end procedure

module procedure new_FullSubspaceBasis_subspace
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  wavevectors = WavevectorBasis( subspace,                  &
                               & frequency,                 &
                               & modes,                     &
                               & qpoints,                   &
                               & maximum_power,             &
                               & potential_expansion_order, &
                               & supercell%sc_size,         &
                               & symmetries                 )
  
  output = FullSubspaceBasis( supercell%sc_size,         &
                            & maximum_power,             &
                            & potential_expansion_order, &
                            & subspace%id,               &
                            & frequency,                 &
                            & wavevectors                )
end procedure

module procedure ground_state_wavefunction
  real(dp) :: mass
  
  real(dp)                  :: coefficient
  type(String), allocatable :: terms(:)
  
  integer :: i,ialloc
  
  ! Calculate the (geometric) average mass.
  mass = product(supercell%atoms%mass())
  mass = mass**(1.0_dp/size(supercell%atoms))
  
  ! Calculate the coefficient.
  coefficient = 1
  allocate(terms(0), stat=ialloc); call err(ialloc)
  do i=1,size(subspace)
    if (subspace%mode_ids(i)==subspace%paired_ids(i)) then
      ! |0_i> = sqrt(sqrt(m*w/pi)) exp(- 1/2 N w (u_i)^2 )
      coefficient = coefficient * (mass*this%frequency/PI)**0.25_dp
      terms = [terms, 'u'//subspace%mode_ids(i)//'^2']
    elseif (subspace%mode_ids(i)<subspace%paired_ids(i)) then
      ! |0_i,0_j> = sqrt(2*m*w/pi) exp(- N w |u_i|^2 )
      coefficient = coefficient * (2*mass*this%frequency/PI)**0.5_dp
      terms = [ terms,                                                    &
              & '2*u'//subspace%mode_ids(i)//'*u'//subspace%paired_ids(i) ]
    endif
  enddo
  
  output = coefficient//'*e^('               // &
     & (-supercell%sc_size*this%frequency/2) // &
     & '*('//join(terms,'+')//'))'
end procedure

module procedure initial_states_FullSubspaceBasis
  associate( gamma_basis =>                                               &
           & this%wavevectors(first(is_int(this%wavevectors%wavevector))) )
    output = gamma_basis%initial_states( subspace,       &
                                       & thermal_energy, &
                                       & anharmonic_data )
  end associate
end procedure

module procedure calculate_states_FullSubspaceBasis
  output = calculate_states( this%wavevectors,        &
                           & subspace,                &
                           & subspace_potential,      &
                           & thermal_energy,          &
                           & state_energy_cutoff,     &
                           & convergence_data,        &
                           & anharmonic_data          )
end procedure

module procedure wavefunction_FullSubspaceBasis
  ! TODO
  output = ''
end procedure

module procedure mode_ids_FullSubspaceBasis
  output = this%wavevectors(1)%mode_ids(subspace,anharmonic_data)
end procedure

module procedure paired_mode_ids_FullSubspaceBasis
  output = this%wavevectors(1)%paired_mode_ids(subspace,anharmonic_data)
end procedure

module procedure inner_product_FullSubspaceBasis
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%inner_product( full_bra,       &
                                            & full_ket,       &
                                            & subspace,       &
                                            & anharmonic_data )
end procedure

module procedure integrate_BasisState_FullSubspaceBasis
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%integrate( full_bra,       &
                                        & monomial,       &
                                        & full_ket,       &
                                        & subspace,       &
                                        & anharmonic_data )
end procedure

module procedure kinetic_energy_FullSubspaceBasis
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_energy( full_bra,       &
                                             & full_ket,       &
                                             & subspace,       &
                                             & anharmonic_data )
end procedure

module procedure harmonic_potential_energy_FullSubspaceBasis
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%harmonic_potential_energy( full_bra,       &
                                                        & full_ket,       &
                                                        & subspace,       &
                                                        & anharmonic_data )
end procedure

module procedure kinetic_stress_FullSubspaceBasis
  type(WavevectorState)              :: full_bra
  type(WavevectorState), allocatable :: full_ket
  
  integer :: i
  
  full_bra = WavevectorState(bra)
  if (present(ket)) then
    full_ket = WavevectorState(ket)
  endif
  
  i = first(this%wavevectors%wavevector == full_bra%wavevector)
  
  output = this%wavevectors(i)%kinetic_stress( full_bra,          &
                                             & full_ket,          &
                                             & subspace,          &
                                             & stress_prefactors, &
                                             & anharmonic_data    )
end procedure

module procedure thermodynamic_data_FullSubspaceBasis
  ! Calculate the thermodynamic properties for the system.
  output = core_shell_thermodynamics( &
              & thermal_energy,       &
              & this%frequency,       &
              & this%supercell_size,  &
              & size(subspace),       &
              & subspace,             &
              & this%wavevectors,     &
              & states,               &
              & subspace_potential,   &
              & subspace_stress,      &
              & stress_prefactors,    &
              & anharmonic_data       )
end procedure

module procedure wavefunctions_FullSubspaceBasis
  type(WavevectorStates), pointer :: full_states
  
  type(String)                    :: ground_state
  type(String), allocatable       :: state_wavefunctions(:)
  type(FullSubspaceWavefunctions) :: wavefunctions
  
  full_states => wavevector_states_pointer(states)
  
  ! Construct the wavefunction of |0>.
  ground_state = this%ground_state_wavefunction(      &
               & subspace,                            &
               & anharmonic_data%anharmonic_supercell )
  
  ! Construct the wavefunctions of each state, in terms of |0>.
  state_wavefunctions = this%wavefunction(  &
     & full_states%states,                  &
     & anharmonic_data%anharmonic_supercell )
  
  ! Pack the wavefunctions into the output,
  !    of concrete type FullSubspaceWavefunctions.
  wavefunctions = FullSubspaceWavefunctions( &
                     & subspace%id,          &
                     & subspace%mode_ids,    &
                     & subspace%paired_ids,  &
                     & ground_state,         &
                     & full_states%energies, &
                     & state_wavefunctions   )
  
  ! Convert the output to abstract class SubspaceWavefunctions.
  output = SubspaceWavefunctionsPointer(wavefunctions)
end procedure

module procedure integrate_BasisStates_FullSubspaceBasis
  integer :: i
  
  output = 0
  do i=1,size(this%wavevectors)
    output = output + this%wavevectors(i)%integrate( states,         &
                                                   & monomial,       &
                                                   & subspace,       &
                                                   & anharmonic_data )
  enddo
end procedure

module procedure free_energy_gradient_FullSubspaceBasis
  output = free_energy_gradient( this%wavevectors,    &
                               & subspace_potential,  &
                               & basis_functions,     &
                               & subspace,            &
                               & states,              &
                               & thermal_energy,      &
                               & state_energy_cutoff, &
                               & anharmonic_data      )
end procedure

module procedure read_FullSubspaceBasis
  integer                            :: supercell_size
  integer                            :: maximum_power
  integer                            :: expansion_order
  integer                            :: subspace_id
  real(dp)                           :: frequency
  type(WavevectorBasis), allocatable :: wavevectors(:)
  
  integer, allocatable :: starting_lines(:)
  integer, allocatable :: ending_lines(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: wavevector_lines(:)
  
  integer :: i,ialloc
  
  select type(this); type is(FullSubspaceBasis)
    line = split_line(input(1))
    supercell_size = int(line(4))
    
    line = split_line(input(2))
    maximum_power = int(line(4))
    
    line = split_line(input(3))
    expansion_order = int(line(4))
    
    line = split_line(input(4))
    subspace_id = int(line(3))
    
    line = split_line(input(5))
    frequency = dble(line(3))
    
    starting_lines = [integer::]
    do i=6,size(input)
      line = split_line(input(i))
      if (size(line)>0) then
        if (line(1)=='Wavevector') then
          starting_lines = [starting_lines, i]
        endif
      endif
    enddo
    
    ending_lines = [starting_lines(2:)-1, size(input)]
    
    allocate(wavevectors(size(starting_lines)), stat=ialloc); call err(ialloc)
    do i=1,size(starting_lines)
      wavevector_lines = [input(:5), input(starting_lines(i):ending_lines(i))]
      wavevectors(i) = WavevectorBasis(wavevector_lines)
    enddo
    
    this = FullSubspaceBasis( supercell_size,  &
                            & maximum_power,   &
                            & expansion_order, &
                            & subspace_id,     &
                            & frequency,       &
                            & wavevectors      )
  class default
    call err()
  end select
end procedure

module procedure write_FullSubspaceBasis
  type(String), allocatable :: wavevector_strings(:)
  
  integer :: i
  
  select type(this); type is(FullSubspaceBasis)
    output = [ 'Supercell size  : '//this%supercell_size,  &
             & 'Maximum power   : '//this%maximum_power,   &
             & 'Expansion order : '//this%expansion_order, &
             & 'Subspace        : '//this%subspace_id,     &
             & 'Frequency       : '//this%frequency        ]
    do i=1,size(this%wavevectors)
      wavevector_strings = str(this%wavevectors(i))
      output = [ output, wavevector_strings(5:) ]
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_FullSubspaceBasis_Strings
  call this%read(input)
end procedure

module procedure new_FullSubspaceBasis_StringArray
  this = FullSubspaceBasis(str(input))
end procedure
end submodule

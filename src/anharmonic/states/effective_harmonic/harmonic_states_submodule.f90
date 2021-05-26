submodule (caesar_harmonic_states_module) caesar_harmonic_states_submodule
  use caesar_effective_harmonic_module
contains

module procedure representation_HarmonicStates
  output = 'harmonic state'
end procedure

module procedure new_HarmonicStates
  this%subspace_id    = subspace_id
  this%frequency      = frequency
  this%thermal_energy = thermal_energy
  this%expectation_cache = ExpectationCache()
end procedure

module procedure new_HarmonicStates_BasisStates
  select type(input); type is(HarmonicStates)
    this = input
  type is(BasisStatesPointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this module procedure
    !    from within this module procedure, so the full name is used instead.
    this = new_HarmonicStates_BasisStates(input%states())
  class default
    call err()
  end select
end procedure

module procedure harmonic_states_pointer
  select type(input); type is(HarmonicStates)
    this => input
  type is(BasisStatesPointer)
    this => harmonic_states_pointer(input%states_pointer())
  class default
    call err()
  end select
end procedure

module procedure read_HarmonicStates
  integer  :: subspace_id
  real(dp) :: frequency
  real(dp) :: thermal_energy
  
  select type(this); type is(HarmonicStates)
    subspace_id = int(token(input(1),3))
    frequency = dble(token(input(2),3))
    thermal_energy = dble(token(input(3),4))
    
    this = HarmonicStates(subspace_id, frequency, thermal_energy)
  class default
    call err()
  end select
end procedure

module procedure write_HarmonicStates
  select type(this); type is(HarmonicStates)
    output = [ 'Subspace       : '//this%subspace_id,   &
             & 'Frequency      : '//this%frequency,     &
             & 'Thermal energy : '//this%thermal_energy ]
  class default
    call err()
  end select
end procedure

module procedure new_HarmonicStates_Strings
  call this%read(input)
end procedure

module procedure new_HarmonicStates_StringArray
  this = HarmonicStates(str(input))
end procedure
end submodule

submodule (caesar_wavevector_state_module) caesar_wavevector_state_submodule
  use caesar_states_module
contains

module procedure startup_wavevector_state
  type(WavevectorState) :: state
  
  call state%startup()
end procedure

module procedure representation_WavevectorState
  output = 'wavevector state'
end procedure

module procedure new_WavevectorState
  this%subspace_id  = subspace_id
  this%wavevector   = wavevector
  this%state_ids    = state_ids
  this%coefficients = coefficients
end procedure

module procedure new_WavevectorState_BasisState
  select type(input); type is(WavevectorState)
    this = input
  type is(BasisStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this module procedure
    !    from within this module procedure, so the full name is used instead.
    this = new_WavevectorState_BasisState(input%state())
  class default
    call err()
  end select
end procedure

module procedure wavevector_state_pointer
  select type(input); type is(WavevectorState)
    this => input
  type is(BasisStatePointer)
    this => wavevector_state_pointer(input%state())
  class default
    call err()
  end select
end procedure

module procedure read_WavevectorState
  integer               :: subspace_id
  type(FractionVector)  :: wavevector
  integer,  allocatable :: state_ids(:)
  real(dp), allocatable :: coefficients(:)
  
  select type(this); type is(WavevectorState)
    subspace_id = int(token(input(1), 3))
    wavevector = FractionVector(join(tokens(input(2),3)))
    state_ids = int(tokens(input(3),3))
    coefficients = dble(tokens(input(4),3))
    
    this = WavevectorState(subspace_id,wavevector,state_ids,coefficients)
  end select
end procedure

module procedure write_WavevectorState
  select type(this); type is(WavevectorState)
    output = [ 'Subspace     : '//this%subspace_id,                      &
             & 'Wavevector   : '//str(this%wavevector),                  &
             & 'States       : '//this%state_ids,                        &
             & 'Coefficients : '//join(this%coefficients, delimiter=' ') ]
  end select
end procedure

module procedure new_WavevectorState_Strings
  call this%read(input)
end procedure

module procedure new_WavevectorState_StringArray
  this = WavevectorState(str(input))
end procedure
end submodule

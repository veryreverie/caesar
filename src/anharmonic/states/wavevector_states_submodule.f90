submodule (caesar_wavevector_states_module) caesar_wavevector_states_submodule
  use caesar_states_module
contains

module procedure new_WavevectorStates
  this%subspace_id = subspace_id
  this%states = states
  this%energies = energies
  if (present(weights)) then
    this%weights = weights
  endif
  this%expectation_cache = ExpectationCache()
end procedure

module procedure new_WavevectorStates_BasisStates
  select type(input); type is(WavevectorStates)
    this = input
  type is(BasisStatesPointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this module procedure
    !    from within this module procedure, so the full name is used instead.
    this = new_WavevectorStates_BasisStates(input%states())
  class default
    call err()
  end select
end procedure

module procedure wavevector_states_pointer
  select type(input); type is(WavevectorStates)
    this => input
  type is(BasisStatesPointer)
    this => wavevector_states_pointer(input%states_pointer())
  class default
    call err()
  end select
end procedure

module procedure representation_WavevectorStates
  output = 'wavevector state'
end procedure

module procedure read_WavevectorStates
  type(StringArray), allocatable :: sections(:)
  
  integer                            :: subspace_id
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  
  integer :: i,ialloc
  
  select type(this); type is(WavevectorStates)
    subspace_id = int(token(input(1),2))
    
    sections = split_into_sections(input(2:), separating_line=repeat('=',50))
    
    allocate( states(size(sections)),   &
            & energies(size(sections)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(sections)
      states(i) = WavevectorState(sections(i)%strings(:size(sections(i))-1))
      
      energies(i) = dble(token(sections(i)%strings(size(sections(i))), 3))
    enddo
    this = WavevectorStates(subspace_id, states, energies)
  class default
    call err()
  end select
end procedure

module procedure write_WavevectorStates
  type(StringArray), allocatable :: sections(:)
  
  integer :: i
  
  select type(this); type is(WavevectorStates)
    sections = [( StringArray([ str(this%states(i)),              &
                &               'Energy : '//this%energies(i) ]), &
                & i=1,                                            &
                & size(this%states)                               )]
    output = [ 'Subspace: '//this%subspace_id,               &
             & str(sections, separating_line=repeat('=',50)) ]
  class default
    call err()
  end select
end procedure

module procedure new_WavevectorStates_Strings
  call this%read(input)
end procedure

module procedure new_WavevectorStates_StringArray
  this = WavevectorStates(str(input))
end procedure
end submodule

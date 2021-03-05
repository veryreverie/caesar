submodule (caesar_subspace_state_module) caesar_subspace_state_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SubspaceStatePointer
  integer :: ialloc
  
  select type(state); type is(SubspaceStatePointer)
    this = state
  class default
    this%representation_ = state%representation()
    allocate( this%state_, source=state, &
            & stat=ialloc); call err(ialloc)
  end select
end procedure

module procedure check_SubspaceStatePointer
  if (.not. allocated(this%state_)) then
    call print_line(CODE_ERROR//': Trying to use a SubspaceStatePointer &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_SubspaceStatePointer
  output = 'pointer'
end procedure

module procedure state_SubspaceStatePointer
  output = this%state_
end procedure

module procedure state_pointer_SubspaceStatePointer
  output => this%state_
end procedure

module procedure mode_ids_SubspaceStatePointer
  call this%check()
  
  output = this%state_%mode_ids()
end procedure

module procedure paired_mode_ids_SubspaceStatePointer
  call this%check()
  
  output = this%state_%paired_mode_ids()
end procedure

module procedure occupation_SubspaceStatePointer
  call this%check()
  
  output = this%state_%occupation()
end procedure

module procedure wavevector_SubspaceStatePointer
  call this%check()
  
  output = this%state_%wavevector(modes,qpoints)
end procedure

module procedure read_SubspaceStatePointer
  type(SubspaceStatePointer), allocatable :: types(:)
  
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(SubspaceStatePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    types = types_SubspaceState()
    i = first([(                                           &
       & types(i)%state_%representation()==representation, &
       & i=1,                                              &
       & size(types)                         )])
    this = types(i)
    
    call this%state_%read(input(2:))
  class default
    call err()
  end select
end procedure

module procedure write_SubspaceStatePointer
  select type(this); type is(SubspaceStatePointer)
    output = [ 'SubspaceState representation: '//this%representation_, &
             & str(this%state_)                                        ]
  end select
end procedure

module procedure new_SubspaceStatePointer_Strings
  call this%read(input)
end procedure

module procedure new_SubspaceStatePointer_StringArray
  this = SubspaceStatePointer(str(input))
end procedure
end submodule

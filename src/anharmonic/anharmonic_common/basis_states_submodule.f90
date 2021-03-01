submodule (caesar_basis_states_module) caesar_basis_states_submodule
  use caesar_anharmonic_common_module
  
  ! An array of all types which extend BasisStates.
  ! This array will be filled in by startup routines.
  type(BasisStatesPointer), allocatable :: TYPES_BasisStates(:)
contains

module procedure startup_BasisStates
  integer :: i
  
  if (.not.allocated(TYPES_BasisStates)) then
    TYPES_BasisStates = [BasisStatesPointer(this)]
  elseif (.not.any([(                                       &
     & this%representation()                                &
     &    == TYPES_BasisStates(i)%states_%representation(), &
     & i=1,                                                 &
     & size(TYPES_BasisStates)                              )])) then
    TYPES_BasisStates = [TYPES_BasisStates, BasisStatesPointer(this)]
  endif
end procedure

module procedure new_BasisStatesPointer
  integer :: ialloc
  
  select type(states); type is(BasisStatesPointer)
    this = states
  class default
    this%representation_ = states%representation()
    allocate( this%states_, source=states, &
            & stat=ialloc); call err(ialloc)
    this%subspace_id = this%states_%subspace_id
  end select
end procedure

module procedure check_BasisStatesPointer
  if (.not. allocated(this%states_)) then
    call print_line(CODE_ERROR//': Trying to use a BasisStatesPointer &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_BasisStatesPointer
  output = 'pointer'
end procedure

module procedure states_BasisStatesPointer
  call this%check()
  
  output = this%states_
end procedure

module procedure states_pointer_BasisStatesPointer
  call this%check()
  
  output => this%states_
end procedure

module procedure read_BasisStatesPointer
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(BasisStatesPointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                        &
       & TYPES_BasisStates(i)%states_%representation()==representation, &
       & i=1,                                                           &
       & size(TYPES_BasisStates)                                        )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_BasisStates(i)%states_%read(input(2:))
    this = BasisStatesPointer(TYPES_BasisStates(i))
  class default
    call err()
  end select
end procedure

module procedure write_BasisStatesPointer
  select type(this); type is(BasisStatesPointer)
    output = [ 'BasisStates representation: '//this%representation_, &
             & str(this%states_)                                     ]
  end select
end procedure

module procedure new_BasisStatesPointer_Strings
  call this%read(input)
end procedure

module procedure new_BasisStatesPointer_StringArray
  this = BasisStatesPointer(str(input))
end procedure
end submodule

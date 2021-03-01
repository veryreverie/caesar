submodule (caesar_basis_state_module) caesar_basis_state_submodule
  use caesar_anharmonic_common_module
  
  ! An array of all types which extend BasisState.
  ! This array will be filled in by startup routines.
  type(BasisStatePointer), allocatable :: TYPES_BasisState(:)
contains

module procedure startup_BasisState
  integer :: i
  
  if (.not.allocated(TYPES_BasisState)) then
    TYPES_BasisState = [BasisStatePointer(this)]
  elseif (.not.any([(                                     &
     & this%representation()                              &
     &    == TYPES_BasisState(i)%state_%representation(), &
     & i=1,                                               &
     & size(TYPES_BasisState)                             )])) then
    TYPES_BasisState = [TYPES_BasisState, BasisStatePointer(this)]
  endif
end procedure

module procedure new_BasisStatePointer
  integer :: ialloc
  
  select type(state); type is(BasisStatePointer)
    this = state
  class default
    this%representation_ = state%representation()
    allocate( this%state_, source=state, &
            & stat=ialloc); call err(ialloc)
    this%subspace_id = this%state_%subspace_id
  end select
end procedure

module procedure check_BasisStatePointer
  if (.not. allocated(this%state_)) then
    call print_line(CODE_ERROR//': Trying to use a BasisStatePointer &
       &before it has been allocated.')
    call err()
  endif
end procedure

module procedure representation_BasisStatePointer
  output = 'pointer'
end procedure

module procedure state_BasisStatePointer
  output = this%state_
end procedure

module procedure state_pointer_BasisStatePointer
  output => this%state_
end procedure

module procedure read_BasisStatePointer
  type(String), allocatable :: line(:)
  
  type(String) :: representation
  
  integer :: i
  
  select type(this); type is(BasisStatePointer)
    line = split_line(input(1))
    representation = line(3)
    
    ! Identify which type corresponds to the representation.
    i = first([(                                                      &
       & TYPES_BasisState(i)%state_%representation()==representation, &
       & i=1,                                                         &
       & size(TYPES_BasisState)                                       )])
    
    ! Read the input into the element of the correct type,
    !    and copy that element into the output.
    call TYPES_BasisState(i)%state_%read(input(2:))
    this = BasisStatePointer(TYPES_BasisState(i))
  class default
    call err()
  end select
end procedure

module procedure write_BasisStatePointer
  select type(this); type is(BasisStatePointer)
    output = [ 'BasisState representation: '//this%representation_, &
             & str(this%state_)                                     ]
  end select
end procedure

module procedure new_BasisStatePointer_Strings
  call this%read(input)
end procedure

module procedure new_BasisStatePointer_StringArray
  this = BasisStatePointer(str(input))
end procedure
end submodule

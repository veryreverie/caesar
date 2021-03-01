submodule (caesar_coupled_states_module) caesar_coupled_states_submodule
  use caesar_states_module
contains

module procedure new_CoupledStates_null
  this%length_ = 0
  this%ids_ = [integer::]
  this%separations_ = [integer::]
end procedure

module procedure new_CoupledStates
  if (size(ids)/=size(separations)) then
    call print_line(CODE_ERROR//': ids and separations do not match.')
    call err()
  endif
  
  this%length_ = size(ids)
  this%ids_ = ids
  this%separations_ = separations
end procedure

module procedure size_CoupledStates
  output = this%length_
end procedure

module procedure id_CoupledStates
  output = this%ids_(index)
end procedure

module procedure ids_CoupledStates
  output = this%ids_(:this%length_)
end procedure

module procedure separation_CoupledStates
  output = this%separations_(index)
end procedure

module procedure separations_CoupledStates
  output = this%separations_(:this%length_)
end procedure

module procedure add_coupling
  integer :: i
  
  if (.not. allocated(this%ids_)) then
    call print_line(CODE_ERROR//': Trying to add a coupling to a &
       &CoupledStates which has not been allocated.')
    call err()
  endif
  
  this%length_ = this%length_ + 1
  if (size(this%ids_)>=this%length_) then
    ! There is space left in the arrays; simply add the coupling.
    this%ids_(this%length_) = id
    this%separations_(this%length_) = separation
  else
    ! There is not sufficient space left in the array.
    ! First double the allocated space, then add the coupling.
    this%ids_ = [this%ids_, id, [(0,i=1,this%length_)]]
    this%separations_ = [this%separations_, separation, [(0,i=1,this%length_)]]
  endif
end procedure

module procedure map_ids
  this%ids_(:this%length_) = id_map(this%ids_(:this%length_))
end procedure

module procedure tensor_product_couplings
  integer :: length
  
  integer :: id
  integer :: separation
  
  integer :: i,j,k,ialloc
  
  length = 0
  do i=1,this%length_
    length = length                                                       &
         & + count( this%separations_(i)+that%separations_(:that%length_) &
         &       <= max_separation                                        &
         &    .and. that%ids_(:that%length_)                              &
         &       <= no_states(unmapped_this%ids_(i))                      )
  enddo
  
  allocate( output%ids_(length),         &
          & output%separations_(length), &
          & stat=ialloc); call err(ialloc)
  output%length_ = length
  
  k = 0
  do i=1,this%length_
    do j=1,that%length_
      if (that%ids_(j)>no_states(unmapped_this%ids_(i))) then
        exit
      endif
      separation = this%separations_(i) + that%separations_(j)
      if (separation>max_separation) then
        cycle
      endif
      id = this%ids_(i) + that%ids_(j)
      k = k+1
      output%ids_(k) = id
      output%separations_(k) = separation
    enddo
  enddo
end procedure

module procedure selected_states_couplings
  integer, allocatable :: state_map(:)
  
  integer :: i,ialloc
  
  ! Generate a list such that state_map(selected_states(i))=i,
  !    and state_map(i)=0 if i is not in selected_states.
  !state_map = [( first(selected_states==i, default=0), &
  !             & i=1,                                  &
  !             & size(input)                           )]
  ! WORKAROUND to avoid memory leak in ifort 19.1.0.166
  allocate(state_map(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    state_map(i) = first(selected_states==i, default=0)
  enddo
  
  ! Generate the couplings.
  allocate(output(size(selected_states)), stat=ialloc); call err(ialloc)
  do i=1,size(selected_states)
    output(i)%i = state_map(input(selected_states(i))%ids())
    output(i)%i = output(i)%i(filter(output(i)%i/=0))
  enddo
end procedure

module procedure read_CoupledStates
  integer, allocatable :: ids(:)
  integer, allocatable :: separations(:)
  
  type(String), allocatable :: line(:)
  type(String), allocatable :: coupling(:)
  
  integer :: i,ialloc
  
  select type(this); type is(CoupledStates)
    line = split_line(input)
    allocate( ids(size(line)),         &
            & separations(size(line)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(line)
      coupling = split_line(line(i),delimiter=':')
      ids(i) = int(coupling(1))
      separations(i) = int(coupling(2))
    enddo
    this = CoupledStates(ids,separations)
  class default
    call err()
  end select
end procedure

module procedure write_CoupledStates
  integer :: i
  
  select type(this); type is(CoupledStates)
    output = join([(this%ids_(i)//':'//this%separations_(i),i=1,size(this))])
  class default
    call err()
  end select
end procedure

module procedure new_CoupledStates_String
  call this%read(input)
end procedure
end submodule

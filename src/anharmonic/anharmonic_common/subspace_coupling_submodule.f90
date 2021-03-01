submodule (caesar_subspace_coupling_module) caesar_subspace_coupling_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SubspaceCoupling
  integer :: ialloc
  
  if (present(ids)) then
    output%ids = ids
  else
    allocate(output%ids(0), stat=ialloc); call err(ialloc)
  endif
end procedure

module procedure concatenate_SubspaceCoupling_DegenerateSubspace
  output%ids = [input%ids,subspace%id]
end procedure

module procedure size_SubspaceCoupling
  output = size(this%ids)
end procedure

module procedure equality_SubspaceCoupling
  integer, allocatable :: this_ids(:)
  integer, allocatable :: that_ids(:)
  
  this_ids = this%ids
  this_ids = this_ids(set(this_ids))
  this_ids = this_ids(sort(this_ids))
  
  that_ids = that%ids
  that_ids = that_ids(set(that_ids))
  that_ids = that_ids(sort(that_ids))
  
  if (size(this_ids)/=size(that_ids)) then
    output = .false.
  else
    output = all(this_ids==that_ids)
  endif
end procedure

module procedure non_equality_SubspaceCoupling
  output = .not. this==that
end procedure

module procedure generate_coupled_subspaces
  type(SubspaceCoupling), allocatable :: temp1(:)
  type(SubspaceCoupling), allocatable :: temp2(:)
  
  integer :: coupling_order
  
  integer :: ialloc
  
  ! Check input.
  if (maximum_coupling_order<1) then
    call print_line(ERROR//': maximum_coupling_order must be at least 1.')
    call quit()
  endif
  
  ! Call the helper module procedure once for each coupling order.
  ! Each call returns an array of results, and these arrays are concatenated
  !    together.
  allocate(output(0), stat=ialloc); call err(ialloc)
  do coupling_order=1,maximum_coupling_order
    ! WORKAROUND: this is done in three steps rather than one to avoid a
    !    compiler bug in ifort 19.0.4.
    temp1 = generate_coupled_subspaces_helper( SubspaceCoupling(), &
                                             & subspaces,          &
                                             & coupling_order      )
    temp2 = [output, temp1]
    output = temp2
  enddo
end procedure

module procedure generate_coupled_subspaces_helper
  type(SubspaceCoupling) :: coupling
  
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(subspaces)
    coupling = coupling_in//subspaces(i)
    if (size(coupling)==coupling_order) then
      output = [output, coupling]
    else
      output = [ output,                                             &
             &   generate_coupled_subspaces_helper( coupling,        &
             &                                      subspaces(i+1:), &
             &                                      coupling_order)  &
             & ]
    endif
  enddo
end procedure

module procedure coupled_subspaces
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = subspaces(first(subspaces%id==this%ids(i)))
  enddo
end procedure

module procedure read_SubspaceCoupling
  select type(this); type is(SubspaceCoupling)
    this = SubspaceCoupling(int(split_line(input)))
  end select
end procedure

module procedure write_SubspaceCoupling
  select type(this); type is(SubspaceCoupling)
    output = join(this%ids)
  end select
end procedure

module procedure new_SubspaceCoupling_String
  call this%read(input)
end procedure
end submodule

submodule (caesar_subspace_coupling_module) caesar_subspace_coupling_submodule
  use caesar_subspaces_module
contains

module procedure new_SubspaceCoupling
  output%ids_ = ids(sort(ids))
  
  ! Check the input.
  if (size(output%ids_)>1) then
    if (any(output%ids_(2:)==output%ids_(:size(output%ids_)-1))) then
      call print_line(ERROR//': A subspace coupling may not contain duplicate &
         &subspaces.')
      call err()
    endif
  endif
end procedure

module procedure ids_SubspaceCoupling
  output = this%ids_
end procedure

module procedure ids_SubspaceCoupling_index
  output = this%ids_(index)
end procedure

module procedure subspaces_SubspaceCoupling
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = subspaces(first(subspaces%id==this%ids_(i)))
  enddo
end procedure

module procedure remove_subspace
  this%ids_ = [this%ids_(:index-1), this%ids_(index+1:)]
end procedure

module procedure read_SubspaceCoupling
  integer, allocatable :: subspace_ids(:)
  
  select type(this); type is(SubspaceCoupling)
    subspace_ids = int(tokens(input, delimiters=['*','s','(',')',' ']))
    this = SubspaceCoupling(subspace_ids)
  class default
    call err()
  end select
end procedure

module procedure write_SubspaceCoupling
  integer :: i
  
  select type(this); type is(SubspaceCoupling)
    output = '('                                                       &
        & // join([('s'//this%ids_(i),i=1,size(this))], delimiter='*') &
        & // ')'
  class default
    call err()
  end select
end procedure

module procedure new_SubspaceCoupling_String
  call this%read(input)
end procedure

module procedure size_SubspaceCoupling
  output = size(this%ids_)
end procedure

module procedure equality_SubspaceCoupling
  if (size(this%ids_)/=size(that%ids_)) then
    output = .false.
  else
    output = all(this%ids_==that%ids_)
  endif
end procedure

module procedure non_equality_SubspaceCoupling
  output = .not. this==that
end procedure

module procedure generate_coupled_subspaces
  type :: IdsAndIndex
    integer, allocatable :: ids(:)
    integer              :: index
  end type
  
  type(IdsAndIndex), allocatable :: old(:)
  type(IdsAndIndex), allocatable :: new(:)
  
  integer :: order        ! The coupling order.
  integer :: no_new_terms ! The number of terms at a given order.
  
  integer :: i,j,k,ialloc
  
  ! Check input.
  if (maximum_coupling_order<1) then
    call print_line(ERROR//': maximum_coupling_order must be at least 1.')
    call quit()
  endif
  
  ! Generate the terms for order=1.
  new = [( IdsAndIndex([subspaces(i)%id],i), i=1, size(subspaces) )]
  ! Calculate the number of terms for order=2.
  no_new_terms = (size(subspaces)*(size(subspaces)-1))/2
  
  do order=2,maximum_coupling_order
    ! Move `new` to `old`, and allocate new to be large enough to hold
    !    `old` plus the terms from the next order.
    old = new
    deallocate(new, stat=ialloc); call err(ialloc)
    allocate(new(size(old)+no_new_terms), stat=ialloc); call err(ialloc)
    
    ! Copy over the terms from the previous order.
    new(:size(old)) = old
    
    ! Generate the terms for the current order,
    !    and calculate how many terms will be generated for the next order.
    no_new_terms = 0
    k = size(old)
    do i=1,size(old)
      do j=old(i)%index+1,size(subspaces)
        k = k+1
        new(k) = IdsAndIndex([old(i)%ids, subspaces(j)%id], j)
        no_new_terms = no_new_terms + (size(subspaces)-j)
      enddo
    enddo
    
    ! Check that the no_new_terms calculation for this order was correct.
    if (k/=size(new)) then
      call print_line(CODE_ERROR//': wrong number of terms calculated.')
      call err()
    endif
  enddo
  
  ! Skip the constructor, as the ids are already in the right format.
  allocate(output(size(new)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%ids_ = new(i)%ids
  enddo
end procedure
end submodule

submodule (caesar_subspace_combination_module) &
   & caesar_subspace_combination_submodule
   use caesar_subspaces_module
contains

module procedure new_SubspaceCombination
  integer, allocatable :: sort_key(:)
  
  if (size(ids)/=size(powers)) then
    call print_line(CODE_ERROR//': IDs and powers do not match.')
    call err()
  elseif (any(powers<1)) then
    call print_line(ERROR//': Every subspace in a combination must have a &
       &power of at least 1.')
    call err()
  endif
  
  sort_key = sort(ids)
  
  this%ids_    = ids(sort_key)
  this%powers_ = powers(sort_key)
  
  if (size(this%ids_)>1) then
    if (any(this%ids_(2:)==this%ids_(:size(this%ids_)-1))) then
      call print_line(ERROR//': A subspace combination may not contain &
         &duplicate subspaces')
      call err()
    endif
  endif
end procedure

module procedure ids_SubspaceCombination
  output = this%ids_
end procedure

module procedure ids_SubspaceCombination_index
  output = this%ids_(index)
end procedure

module procedure powers_SubspaceCombination
  output = this%powers_
end procedure

module procedure powers_SubspaceCombination_index
  output = this%powers_(index)
end procedure

module procedure subspaces_SubspaceCombination
  integer :: i
  
  output = [( subspaces(first(subspaces%id==this%ids_(i))), i=1, size(this) )]
end procedure

module procedure is_subsidiary_of
  integer :: i,j
  
  j = 0
  do i=1,size(this)
    j = j + first_equivalent( that%ids_(j+1:), &
                            & this%ids_(i),    &
                            & default=-j       )
    if (j==0) then
      output = .false.
      return
    elseif (that%powers_(j)<this%powers(i)) then
      output = .false.
      return
    endif
  enddo
  
  output = .true.
end procedure

module procedure read_SubspaceCombination
  integer, allocatable :: ids(:)
  integer, allocatable :: powers(:)
  
  type(String), allocatable :: subspace_strings(:)
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceCombination)
    subspace_strings = tokens(input, delimiters=['*','s','(',')'])
    allocate( ids(size(subspace_strings)),    &
            & powers(size(subspace_strings)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(subspace_strings)
      line = tokens(subspace_strings(i), delimiter='^')
      ids(i) = int(line(1))
      powers(i) = int(line(2))
    enddo
    this = SubspaceCombination(ids, powers)
  end select
end procedure

module procedure write_SubspaceCombination
  type(String), allocatable :: subspace_strings(:)
  
  integer :: i
  
  select type(this); type is(SubspaceCombination)
    subspace_strings = [( 's'//this%ids_(i)//'^'//this%powers_(i), &
                        & i=1,                                     &
                        & size(this)                               )]
    output = '('//join(subspace_strings, delimiter='*')//')'
  end select
end procedure

module procedure new_SubspaceCombination_String
  call this%read(input)
end procedure

module procedure size_SubspaceCombination
  output = size(this%ids_)
end procedure

module procedure equality_SubspaceCombination_SubspaceCombination
  if (size(this)/=size(that)) then
    output = .false.
  else
    output = all(this%ids_==that%ids_) .and. all(this%powers_==that%powers_)
  endif
end procedure

module procedure non_equality_SubspaceCombination_SubspaceCombination
  output = .not. this==that
end procedure

module procedure lt_SubspaceCombination_SubspaceCombination
  integer :: i
  
  if (sum(this%powers_)<sum(that%powers_)) then
    output = .true.
  elseif (sum(this%powers_)>sum(that%powers_)) then
    output = .false.
  else
    do i=1,min(size(this%ids_),size(that%ids_))
      if (this%ids_(i)<that%ids_(i)) then
        output = .true.
        return
      elseif (this%ids_(i)>that%ids_(i)) then
        output = .false.
        return
      elseif (this%powers_(i)<that%powers_(i)) then
        output = .false.
        return
      elseif (this%powers_(i)>that%powers_(i)) then
        output = .true.
        return
      endif
    enddo
    
    output = .false.
  endif
end procedure

module procedure le_SubspaceCombination_SubspaceCombination
  output = this==that .or. this<that
end procedure

module procedure gt_SubspaceCombination_SubspaceCombination
  output = .not. this<=that
end procedure

module procedure ge_SubspaceCombination_SubspaceCombination
  output = .not. this<that
end procedure

module procedure generate_subspace_combinations
  type :: SubspacePowers
    integer, allocatable :: powers(:)
    integer              :: total_power
    integer              :: next_starting_power
    integer              :: next_stopping_power
  end type
  
  type(SubspacePowers), allocatable :: old(:)
  type(SubspacePowers), allocatable :: new(:)
  
  integer :: no_new_terms
  
  integer :: i,j,k,l,ialloc
  
  ! Handle the case of no subspace couplings separately.
  if (size(subspace_coupling)==0 .and. minimum_power>0) then
    output = [SubspaceCombination::]
    return
  endif
  
  ! Seed the calculation with powers=0.
  if (size(subspace_coupling)==1) then
    new = [SubspacePowers(                                            &
           & powers              = [(0,i=1,size(subspace_coupling))], &
           & total_power         = 0,                                 &
           & next_starting_power = max(1,minimum_power),              &
           & next_stopping_power = maximum_power                      )]
  else
    new = [SubspacePowers(                                                 &
           & powers              = [(0,i=1,size(subspace_coupling))],      &
           & total_power         = 0,                                      &
           & next_starting_power = 1,                                      &
           & next_stopping_power = maximum_power-size(subspace_coupling)+1 )]
  endif
  no_new_terms = new(1)%next_stopping_power-new(1)%next_starting_power+1
  
  ! Loop over the subspaces.
  ! For each subspace, loop over all combinations from the previous subspace,
  !    and add a combination for each allowed power of the new subspace.
  do i=1,size(subspace_coupling)
    ! Move `new` to `old`, and allocate new to be large enough to hold
    !    the terms from the next order.
    old = new
    deallocate(new, stat=ialloc); call err(ialloc)
    allocate(new(no_new_terms), stat=ialloc); call err(ialloc)
    
    ! Generate the terms for the current subspace,
    !    and calculate how many terms will be generated for the next subspace.
    l = 0
    no_new_terms = 0
    do j=1,size(old)
      do k=old(j)%next_starting_power,old(j)%next_stopping_power
        l = l+1
        new(l) = old(j)
        new(l)%powers(i) = k
        new(l)%total_power = new(l)%total_power+k
        if (i+1<size(subspace_coupling)) then
          new(l)%next_starting_power = 1
          new(l)%next_stopping_power = maximum_power-new(l)%total_power
        elseif (i+1==size(subspace_coupling)) then
          new(l)%next_starting_power = max(1,minimum_power-new(l)%total_power)
          new(l)%next_stopping_power = maximum_power-new(l)%total_power &
                                   & - size(subspace_coupling)+i+1
        else
          ! next_*_power will not be used, so does not need calculating.
          continue
        endif
        
        no_new_terms = no_new_terms &
                   & + new(l)%next_stopping_power-new(l)%next_starting_power+1
      enddo
    enddo
    
    ! Check that the no_new_terms calculation for this subspace was correct.
    if (l/=size(new)) then
      call print_line(CODE_ERROR//': wrong number of terms calculated.')
      call err()
    endif
  enddo
  
  ! Skip the constructor,
  !    as the ids and powers are already in the right format.
  allocate(output(size(new)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%ids_    = subspace_coupling%ids()
    output(i)%powers_ = new(i)%powers
  enddo
end procedure
end submodule

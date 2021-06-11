submodule (caesar_combination_qpoint_combination_module) &
   & caesar_combination_qpoint_combination_submodule
  use caesar_stars_module
contains

module procedure new_CombinationQpointCombination
  if (size(subspace_combination)/=size(qpoint_combinations)) then
    call print_line(ERROR//': subspace_combination and qpoint_combinations do &
       &not match.')
    call err()
  elseif (any( subspace_combination%powers()     &
          & /= qpoint_combinations%total_power() )) then
    call print_line(ERROR//': subspace_combination and qpoint_combinations do &
       &not match.')
    call print_line(subspace_combination)
    call print_lines(qpoint_combinations)
    call err()
  endif
  
  this%subspace_combination = subspace_combination
  this%qpoint_combinations  = qpoint_combinations
end procedure

module procedure total_power_CombinationQpointCombination
  output = sum(this%qpoint_combinations%total_power())
end procedure

module procedure wavevector_CombinationQpointCombination
  integer :: i
  
  if (size(this%qpoint_combinations)==0) then
    output = fracvec(zeroes(3))
  else
    output = sum([( this%qpoint_combinations(i)%wavevector(qpoints), &
                  & i=1,                                             &
                  & size(this%qpoint_combinations)                   )])
  endif
end procedure

module procedure complex_monomials_CombinationQpointCombination
  type(ComplexMonomial), allocatable :: old(:)
  type(ComplexMonomial), allocatable :: new(:)
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  type(ComplexMonomial), allocatable :: subspace_monomials(:)
  
  integer :: i,j,k,l,ialloc
  
  new = [ComplexMonomial((1.0_dp,0.0_dp), [ComplexUnivariate::])]
  do i=1,size(this%qpoint_combinations)
    ! Get the modes in the i'th subspace.
    subspace_modes = modes(                                          &
       & filter(modes%subspace_id==this%subspace_combination%ids(i)) )
    
    ! Calcualte the monomials corresponding to qpoint_combinations(i)
    !    in the i'th subspace.
    subspace_monomials = this%qpoint_combinations(i)%complex_monomials( &
                                                       & subspace_modes )
    
    ! Copy the array "new" to "old", and re-allocate "new" to be large enough
    !    to store the monomials for this iteration.
    old = new
    deallocate(new, stat=ialloc); call err(ialloc)
    allocate( new(size(old)*(size(subspace_monomials))), &
            & stat=ialloc); call err(ialloc)
    
    ! Loop over "old" and "qpoint_monomials",
    !    and construct all monomials which are products of the two.
    l = 0
    do j=1,size(old)
      do k=1,size(subspace_monomials)
        l = l+1
        new(l) = old(j)*subspace_monomials(k)
      enddo
    enddo
  enddo
  
  output = new
end procedure

module procedure read_CombinationQpointCombination
  type(SubspaceCombination)            :: subspace_combination
  type(QpointCombination), allocatable :: qpoint_combinations(:)
  
  type(String), allocatable :: line(:)
  select type(this); type is(CombinationQpointCombination)
    subspace_combination = SubspaceCombination(token(input(1), 6))
    line = tokens(input(2))
    if (size(subspace_combination)==0) then
      ! If size(subspace_combination)==0 then the second line should be '()'.
      if (size(line)/=1 .or. any(line/=str('()'))) then
        call print_line(ERROR//': Unable to parse &
           &CombinationQpointCombination.')
        call print_lines(input)
        call err()
      else
        qpoint_combinations = [QpointCombination::]
      endif
    else
      ! If size(subspace_combination)/=0 then the second line should be the
      !    q-point combinations, separated by ' * '.
      ! This is parsed by separating the line into tokens by spaces,
      !    and ignoring the '*' tokens.
      if (size(line)/=size(subspace_combination)*2-1) then
        call print_line(ERROR//': Unable to parse &
           &CombinationQpointCombination.')
        call print_lines(input)
        call err()
      elseif (any(line(2::2)/=str('*'))) then
        call print_line(ERROR//': Unable to parse &
           &CombinationQpointCombination.')
        call print_lines(input)
        call err()
      else
        qpoint_combinations = QpointCombination(line(1::2))
      endif
    endif
    this = CombinationQpointCombination( subspace_combination, &
                                       & qpoint_combinations   )
  class default
    call err()
  end select
end procedure

module procedure write_CombinationQpointCombination
  select type(this); type is(CombinationQpointCombination)
    output = [ 'q-point combination in subspace combination '// &
                              & this%subspace_combination//' :' ]
    if (size(this%qpoint_combinations)==0) then
      output = [output, str('()')]
    else
      output = [output, join(this%qpoint_combinations, ' * ')]
    endif
  class default
    call err()
  end select
end procedure

module procedure new_CombinationQpointCombination_Strings
  call this%read(input)
end procedure

module procedure new_CombinationQpointCombination_StringArray
  call this%read(str(input))
end procedure

module procedure eq_CombinationQpointCombination_CombinationQpointCombination
  if (size(this%qpoint_combinations)==size(that%qpoint_combinations)) then
    output = this%subspace_combination==that%subspace_combination &
     & .and. all(this%qpoint_combinations==that%qpoint_combinations)
  else
    output = .false.
  endif
end procedure

module procedure ne_CombinationQpointCombination_CombinationQpointCombination
  output = .not. this==that
end procedure

module procedure lt_CombinationQpointCombination_CombinationQpointCombination
  integer :: i
  
  if (this%subspace_combination<that%subspace_combination) then
    output = .true.
  elseif (this%subspace_combination>that%subspace_combination) then
    output = .false.
  else
    do i=1,size(this%qpoint_combinations)
      if (this%qpoint_combinations(i)<that%qpoint_combinations(i)) then
        output = .true.
        return
      elseif (this%qpoint_combinations(i)>that%qpoint_combinations(i)) then
        output = .false.
        return
      endif
    enddo
    
    output = .false.
  endif
end procedure

module procedure le_CombinationQpointCombination_CombinationQpointCombination
  output = this==that .or. this<that
end procedure

module procedure gt_CombinationQpointCombination_CombinationQpointCombination
  output = .not. this<=that
end procedure

module procedure ge_CombinationQpointCombination_CombinationQpointCombination
  output = .not. this<that
end procedure

module procedure conjg_CombinationQpointCombination
  output%subspace_combination = this%subspace_combination
  output%qpoint_combinations = conjg(this%qpoint_combinations)
end procedure

module procedure operate_Group_CombinationQpointCombination
  output = CombinationQpointCombination(     &
     & this%subspace_combination,            &
     & qpoint_group*this%qpoint_combinations )
end procedure

module procedure generate_combination_qpoint_combinations
  type :: CombinationData
    type(CombinationQpointCombination) :: combination
    type(FractionVector)               :: wavevector
    integer, allocatable               :: qpoint_ids(:)
  end type
  
  type(CombinationData), allocatable :: old(:)
  type(CombinationData), allocatable :: new(:)
  
  type(FractionVector), allocatable :: wavevectors(:)
  
  type(QpointCombination) :: combination
  
  integer :: i,j,k,l,ialloc
  
  associate(qpoint_stars => qpoint_star_product%qpoint_stars)
    ! If there are no subspaces, or any subspace contains no stars,
    !    then the output is an empty array.
    if (size(qpoint_stars)==0) then
      output = [CombinationQpointCombination::]
      return
    elseif (any([(size(qpoint_stars(i)),i=1,size(qpoint_stars))]==0)) then
      output = [CombinationQpointCombination::]
      return
    endif
    
    ! Seed the loop with a CombinationQpointCombination whose
    !    qpoint_combinations are not filled in.
    allocate(new(1), stat=ialloc); call err(ialloc)
    new(1)%combination%subspace_combination = &
       & qpoint_star_product%subspace_combination
    allocate( new(1)%combination%qpoint_combinations(size(qpoint_stars)), &
            & stat=ialloc); call err(ialloc)
    new(1)%wavevector = fracvec(zeroes(3))
    new(1)%qpoint_ids = [integer::]
    
    ! Loop over the subspaces.
    ! For each subspace, loop over all previous q-point combinations,
    !    and combine each with each combination in the new subspace.
    do i=1,size(qpoint_stars)
      ! Copy `new` to `old`, and re-allocate new to be large enough to hold
      !    the new combinations.
      old = new
      deallocate(new, stat=ialloc); call err(ialloc)
      allocate( new(size(old)*size(qpoint_stars(i))), &
              & stat=ialloc); call err(ialloc)
      
      ! Calculate the wavevectors for the q-point combinations in the i'th
      !    subspace.
      if (set_default(conserve_momentum,.false.)) then
        wavevectors = qpoint_stars(i)%wavevectors(qpoints)
      endif
      
      ! Loop over the q-point combinations for the first i-1 subspaces,
      !    and the q-point combinations for the i'th subspace,
      !    and generate the products of the two.
      l = 0
      do j=1,size(old)
        do k=1,size(qpoint_stars(i))
          l = l+1
          new(l) = old(j)
          combination = qpoint_stars(i)%combinations(k)
          new(l)%combination%qpoint_combinations(i) = combination
          if (set_default(conserve_momentum,.false.)) then
            new(l)%wavevector = new(l)%wavevector + wavevectors(k)
            if (       i==size(qpoint_stars)              &
               & .and. .not. is_int(new(l)%wavevector) ) then
              l = l-1
              cycle
            endif
          endif
          if (present(max_qpoint_coupling)) then
            new(l)%qpoint_ids = merge_qpoint_ids( old(j)%qpoint_ids, &
                                                & combination        )
            if (size(new(l)%qpoint_ids)>max_qpoint_coupling) then
              l = l-1
              cycle
            endif
          endif
        enddo
      enddo
      
      output = new(:l)%combination
    enddo
  end associate
contains
  ! Extract the q-point ids from `combination`, and merge them with `old_ids`
  !    to give a combined list of q-point ids in ascending order.
  ! e.g.
  !    merge([3,5,7], QpointCombination('(q1^4)*(q5^2)*(q8^4)') = [1,3,5,7,8]
  function merge_qpoint_ids(old_ids,combination) result(output)
    integer,                 intent(in)  :: old_ids(:)
    type(QpointCombination), intent(in)  :: combination
    integer, allocatable                 :: output(:)
    
    type(QpointPower), allocatable :: qpoint_powers(:)
    integer,           allocatable :: new_ids(:)
    
    integer :: i,j,k,ialloc
    
    qpoint_powers = combination%qpoints()
    new_ids = qpoint_powers%id()
    
    allocate( output(size(old_ids)+size(new_ids)), &
            & stat=ialloc); call err(ialloc)
    i = 1
    j = 1
    k = 0
    do while(i<=size(old_ids) .or. j<=size(new_ids))
      k = k+1
      if (i>size(old_ids)) then
        output(k) = new_ids(j)
        j = j+1
      elseif (j>size(new_ids)) then
        output(k) = old_ids(i)
        i = i+1
      elseif (old_ids(i)<new_ids(j)) then
        output(k) = old_ids(i)
        i = i+1
      elseif (new_ids(j)<old_ids(i)) then
        output(k) = new_ids(j)
        j = j+1
      else
        output(k) = old_ids(i)
        i = i+1
        j = j+1
      endif
    enddo
    output = output(:k)
  end function
end procedure
end submodule

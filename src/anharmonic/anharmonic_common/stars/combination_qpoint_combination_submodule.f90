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
  type :: CombinationAndWavevector
    type(CombinationQpointCombination) :: combination
    type(FractionVector)               :: wavevector
  end type
  
  type(CombinationAndWavevector), allocatable :: old(:)
  type(CombinationAndWavevector), allocatable :: new(:)
  
  type(FractionVector), allocatable :: wavevectors(:)
  
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
          new(l)%combination%qpoint_combinations(i) = &
             & qpoint_stars(i)%combinations(k)
          if (set_default(conserve_momentum,.false.)) then
            new(l)%wavevector = new(l)%wavevector + wavevectors(k)
            if (       i==size(qpoint_stars)              &
               & .and. .not. is_int(new(l)%wavevector) ) then
              l = l-1
            endif
          endif
        enddo
      enddo
    enddo
    
    output = new(:l)%combination
  end associate
end procedure
end submodule

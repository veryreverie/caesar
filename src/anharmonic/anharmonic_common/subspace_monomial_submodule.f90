submodule (caesar_subspace_monomial_module) caesar_subspace_monomial_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_SubspaceMonomial
  integer :: ialloc
  
  allocate( this%ids(0),    &
          & this%powers(0), &
          & stat=ialloc); call err(ialloc)
end procedure

module procedure new_SubspaceMonomial_ids_powers
  if (size(ids)/=size(powers)) then
    call print_line(CODE_ERROR//': IDs and powers do not match.')
    call err()
  endif
  
  this%ids = ids
  this%powers = powers
end procedure

module procedure new_SubspaceMonomial_DegenerateSubspaces
  integer :: i
  
  this%ids = subspaces%id
  this%ids = this%ids(sort(this%ids))
  this%ids = this%ids(set(this%ids))
  this%powers = [(count(subspaces%id==this%ids(i)), i=1, size(subspaces))]
end procedure

module procedure concatenate_SubspaceMonomial_DegenerateSubspace
  integer :: i
  
  output = this
  
  i = first(this%ids==subspace%id, default=0)
  if (i==0) then
    output%ids = [output%ids, subspace%id]
    output%powers = [output%powers, 1]
  else
    output%powers(i) = output%powers(i) + 1
  endif
end procedure

module procedure size_SubspaceMonomial
  output = size(this%ids)
end procedure

module procedure equality_SubspaceMonomial_SubspaceMonomial
  if (size(this)/=size(that)) then
    output = .false.
  else
    output = all(this%ids==that%ids) .and. all(this%powers==that%powers)
  endif
end procedure

module procedure non_equality_SubspaceMonomial_SubspaceMonomial
  output = .not. this==that
end procedure

module procedure coupled_subspaces
  integer :: i
  
  output = [( subspaces(first(subspaces%id==this%ids(i))), i=1, size(this) )]
end procedure

module procedure is_subsidiary_of
  integer :: i,j
  
  ! Loop over the ids in this, checking if each is in that.
  do i=1,size(this)
    j = first(that%ids==this%ids(i), default=0)
    if (j==0) then
      output = .false.
      return
    elseif (this%powers(i)>that%powers(j)) then
      output = .false.
      return
    endif
  enddo
  
  output = .true.
end procedure

module procedure generate_subspace_monomials
  type(DegenerateSubspace), allocatable :: coupled_subspaces(:)
  
  ! Check inputs.
  if (size(subspace_coupling)==0) then
    call print_line(CODE_ERROR//': Empty subspace coupling.')
    call err()
  elseif (minimum_expansion_order<0) then
    call print_line(ERROR//': minimum_expansion_order must be non-negative.')
    call quit()
  elseif (   maximum_expansion_order                               &
         & < min(minimum_expansion_order, size(subspace_coupling)) ) then
    call print_line(ERROR//': maximum_expansion_order must be at least as &
       &large as minimum_expansion_order and the size of the coupling.')
    call quit()
  endif
  
  ! Retrieve coupled subspaces from subspace coupling.
  coupled_subspaces = subspace_coupling%coupled_subspaces(subspaces)
  
  ! Call helper module procedure to generate monomials.
  output = generate_subspace_monomials_helper( coupled_subspaces,       &
                                             & minimum_expansion_order, &
                                             & maximum_expansion_order  )
end procedure

module procedure generate_subspace_monomials_helper
  type(DegenerateSubspace)              :: first_subspace
  type(DegenerateSubspace), allocatable :: remaining_subspaces(:)
  
  type(SubspaceMonomial) :: monomial
  
  integer :: ialloc
  
  ! If there are no more subspaces to append, return the monomial,
  !    but only if its size is at least minimum_expansion_order.
  ! Size 1 monomials are ignored because they correspond to linear terms in the
  !    potential, which are zero because the structure is geometry optimised.
  if (size(coupled_subspaces)==0) then
    if (.not. present(monomial_in)) then
      call print_line(CODE_ERROR//': Empty subspace coupling.')
      call err()
    elseif (sum(monomial_in%powers)<minimum_expansion_order) then
      allocate(output(0), stat=ialloc); call err(ialloc)
    else
      output = [monomial_in]
    endif
    
    return
  endif
  
  ! Split the coupled subspaces into the first subspace, which will be handled
  !    by this call of the module procedure, and all the rest, which will be handled
  !    by a call.
  first_subspace      = coupled_subspaces(1)
  remaining_subspaces = coupled_subspaces(2:)
  
  ! Append copies of the first subspace as many times as still leaves space
  !    for at least one copy of every remaining subspace, and then call this
  !    module procedure to handle the next subspace along.
  if (present(monomial_in)) then
    monomial = monomial_in
  else
    monomial = SubspaceMonomial()
  endif
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do while(    sum(monomial%powers)+size(coupled_subspaces) &
          & <= maximum_expansion_order                      )
    monomial = monomial // first_subspace
    output = [ output,                                                      &
           &   generate_subspace_monomials_helper( remaining_subspaces,     &
           &                                       minimum_expansion_order, &
           &                                       maximum_expansion_order, &
           &                                       monomial)                &
           & ]
  enddo
end procedure

module procedure generate_complex_monomials
  type(ComplexUnivariate) :: zero_univariate(0)
  type(ComplexMonomial)   :: root
  
  type(DegenerateSubspace)       :: subspace
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  type(ComplexMonomial), allocatable :: subspace_monomials(:)
  
  integer :: i,j,k,ialloc
  
  if (size(this)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  root = ComplexMonomial( coefficient = cmplx(1.0_dp,0.0_dp,dp), &
                        & modes       = zero_univariate          )
  
  do i=1,size(this)
    subspace = subspaces(first(subspaces%id==this%ids(i)))
    subspace_modes = subspace%modes(modes)
    subspace_monomials = generate_subspace_complex_monomials( subspace_modes, &
                                                            & this%powers(i), &
                                                            & root            )
    if (conserve_subspace_momentum) then
      subspace_monomials = subspace_monomials(filter( subspace_monomials, &
                                                    & conserves_momentum  ))
    endif
    
    if (i==1) then
      output = subspace_monomials
    else
      output = [(                                                            &
         & (output(k)*subspace_monomials(j), j=1, size(subspace_monomials)), &
         & k=1,                                                              &
         & size(output)                                                      )]
    endif
  enddo
  
  if (conserve_momentum) then
    output = output(filter(output, conserves_momentum))
  endif
contains
  ! Lambda for checking if a monomial conserves momentum.
  ! Captures:
  !    - modes
  !    - qpoints
  ! N.B. a monomial given by
  !    prod_i (u_{q_i,i})^{n_i}
  ! conserves momentum iff sum_i n_i q_i is a G-vector.
  function conserves_momentum(input) result(output)
    class(*), intent(in) :: input
    logical              :: output
    select type(input); type is(ComplexMonomial)
      output = input%wavevector(modes,qpoints)==fracvec(zeroes(3))
    end select
  end function
end procedure

module procedure generate_paired_monomials
  type(ComplexMonomial), allocatable :: complex_monomials(:)
  
  complex_monomials = generate_complex_monomials( this,                      &
                                                & maximum_coupling_order,    &
                                                & subspaces,                 &
                                                & modes,                     &
                                                & qpoints,                   &
                                                & conserve_momentum,         &
                                                & conserve_subspace_momentum )
  
  output = PairedMonomial(complex_monomials)
end procedure

module procedure generate_subspace_complex_monomials
  integer :: i
  
  if (size(modes)==0) then
    if (power>0) then
      call print_line(CODE_ERROR//': Monomial requires further modes, but no &
         &modes remain to append.')
      call err()
    endif
  endif
  
  if (power==0) then
    output = [root]
    output%coefficient = sqrt(no_permutations(output(1)))
  elseif (size(modes)==1) then
    output = [root*ComplexUnivariate(mode=modes(1), power=power)]
    output%coefficient = sqrt(no_permutations(output(1)))
  else
    output = [generate_subspace_complex_monomials(modes(2:),power,root)]
    do i=1,power
      output = [ output,                                            &
               & generate_subspace_complex_monomials(               &
               &   modes(2:),                                       &
               &   power-i,                                         &
               &   root*ComplexUnivariate(mode=modes(1), power=i) ) ]
    enddo
  endif
end procedure

module procedure no_permutations
  integer, allocatable :: powers(:)
  
  integer :: i,ialloc
  
  allocate(powers(0), stat=ialloc); call err(ialloc)
  do i=1,size(input)
    powers = [powers, input%power(i)]
    if (input%id(i)/=input%paired_id(i)) then
      powers = [powers, input%paired_power(i)]
    endif
  enddo
  
  output = real_multinomial(sum(powers), powers)
end procedure

module procedure read_SubspaceMonomial
  integer, allocatable :: ids(:)
  integer, allocatable :: powers(:)
  
  type(String), allocatable :: subspace_strings(:)
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(SubspaceMonomial)
    subspace_strings = split_line(input, delimiter='*')
    allocate( ids(size(subspace_strings)),    &
            & powers(size(subspace_strings)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(subspace_strings)
      line = split_line(subspace_strings(i), delimiter='^')
      ids(i) = int(slice(line(1),2,len(line(1))))
      powers(i) = int(slice(line(2),1,len(line(2))-1))
    enddo
    this = SubspaceMonomial(ids, powers)
  end select
end procedure

module procedure write_SubspaceMonomial
  type(String), allocatable :: subspace_strings(:)
  
  integer :: i
  
  select type(this); type is(SubspaceMonomial)
    subspace_strings = [( '(s'//this%ids(i)//'^'//this%powers(i)//')', &
                        & i=1,                                         &
                        & size(this)                                   )]
    output = join(subspace_strings, delimiter='*')
  end select
end procedure

module procedure new_SubspaceMonomial_String
  call this%read(input)
end procedure
end submodule

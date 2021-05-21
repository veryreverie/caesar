submodule (caesar_combination_qpoint_star_module) &
   & caesar_combination_qpoint_star_submodule
  use caesar_stars_module
contains

module procedure new_CombinationQpointStar
  if (size(qpoint_combinations)==0) then
    call print_line(ERROR//': a q-point star must contain at least one &
       &combination.')
    call err()
  elseif (any( qpoint_combinations%subspace_combination    &
      & /= qpoint_combinations(1)%subspace_combination )) then
    call print_line(ERROR//': q-point combinations correspond to different &
       &subspace combinations.')
    call err()
  endif
  
  this%subspace_combination = qpoint_combinations(1)%subspace_combination
  this%qpoint_combinations  = qpoint_combinations(sort( qpoint_combinations, &
                                                      & compare_combinations ))
contains
  function compare_combinations(lhs,rhs) result(output)
    class(*), intent(in) :: lhs
    class(*), intent(in) :: rhs
    logical              :: output
    
    select type(lhs); type is(CombinationQpointCombination)
      select type(rhs); type is(CombinationQpointCombination)
        output = lhs<rhs
      end select
    end select
  end function
end procedure

module procedure complex_monomials_CombinationQpointStar
  integer :: i
  
  output = [( this%qpoint_combinations(i)%complex_monomials(modes), &
            & i=1,                                                  &
            & size(this%qpoint_combinations)                        )]
end procedure

module procedure read_CombinationQpointStar
  type(SubspaceCombination)                       :: subspace_combination
  type(CombinationQpointCombination), allocatable :: qpoint_combinations(:)
  
  integer :: i,ialloc
  
  select type(this); type is(CombinationQpointStar)
    subspace_combination = SubspaceCombination(token(input(1), 6))
    allocate( qpoint_combinations(size(input)-1), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(input)-1
      qpoint_combinations(i) = CombinationQpointCombination(input([1,i+1]))
    enddo
    this = CombinationQpointStar(qpoint_combinations)
  class default
    call err()
  end select
end procedure

module procedure write_CombinationQpointStar
  type(String), allocatable :: lines(:)
  
  integer :: i
  select type(this); type is(CombinationQpointStar)
    output = [ 'q-point star in subspace combination '// &
             &    this%subspace_combination//' :'        ]
    do i=1,size(this%qpoint_combinations)
      lines = str(this%qpoint_combinations(i))
      output = [output, lines(2)]
    enddo
  class default
    call err()
  end select
end procedure

module procedure new_CombinationQpointStar_Strings
  call this%read(input)
end procedure

module procedure new_CombinationQpointStar_StringArray
  call this%read(str(input))
end procedure

module procedure generate_combination_qpoint_stars
  type(CombinationQpointCombination), allocatable :: combinations(:)
  
  integer              :: no_stars
  integer, allocatable :: star_sizes(:)
  integer, allocatable :: star_id(:)
  integer, allocatable :: id_in_star(:)
  
  type(CombinationQpointCombination) :: combination
  
  integer :: i,j,k,ialloc
  
  if (       set_default(conserve_momentum, .false.) &
     & .and. .not. present(qpoints)                  ) then
    call print_line(CODE_ERROR//': If conserve_momentum is true then qpoints &
       &must be present.')
    call err()
  endif
  
  combinations = generate_combination_qpoint_combinations( &
                                    & qpoint_star_product, &
                                    & conserve_momentum,   &
                                    & qpoints              )
  
  if (size(combinations)==0) then
    output = [CombinationQpointStar::]
    return
  endif
  
  allocate( star_sizes(size(combinations)), &
          & star_id(size(combinations)),    &
          & id_in_star(size(combinations)), &
          & stat=ialloc); call err(ialloc)
  star_id = 0
  
  ! Process the combinations to identify the q-point stars.
  !  - There are `no_stars` stars.
  !  - The star with id `i` contains `star_sizes(i)` combinations.
  !  - Combination `combinations(i)` is the `id_in_star(i)` combination
  !       in the star with id `star_id(i)`.
  no_stars = 1
  star_sizes(1) = 1
  star_id(1) = 1
  id_in_star(1) = 1
  
  ! Process each combination (after the first) in order.
  do i=2,size(combinations)
    do j=0,size(qpoint_groups)
      ! Construct a combination which is equivalent to `combinations(i)`,
      !    either by conjugation or by a q-point group.
      if (j==0) then
        combination = conjg(combinations(i))
      else
        combination = qpoint_groups(j) * combinations(i)
      endif
      
      ! If `combinations(i)` is in an existing star, find that star and add
      !    the combination to it.
      if (combination<combinations(i)) then
        k = first(combination==combinations(:i))
        star_sizes(star_id(k)) = star_sizes(star_id(k)) + 1
        star_id(i) = star_id(k)
        id_in_star(i) = star_sizes(star_id(i))
        exit
      endif
    enddo
    
    ! If `combinations(i)` is not in an existing star, create a new star.
    if (star_id(i)==0) then
      no_stars = no_stars + 1
      star_sizes(no_stars) = 1
      star_id(i) = no_stars
      id_in_star(i) = 1
    endif
  enddo
  
  ! Construct the actual q-point stars from the processed q-point combinations.
  allocate(output(no_stars), stat=ialloc); call err(ialloc)
  do i=1,no_stars
    output(i)%subspace_combination = qpoint_star_product%subspace_combination
    allocate( output(i)%qpoint_combinations(star_sizes(i)), &
            & stat=ialloc); call err(ialloc)
  enddo
  do i=1,size(combinations)
    output(star_id(i))%qpoint_combinations(id_in_star(i)) = combinations(i)
  enddo
end procedure
end submodule

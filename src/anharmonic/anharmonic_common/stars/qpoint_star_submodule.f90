submodule (caesar_qpoint_star_module) caesar_qpoint_star_submodule
  use caesar_stars_module
contains

module procedure new_QpointStar
  if (size(combinations)==0) then
    call print_line(CODE_ERROR//': A QpointStar must contain at least one &
       &QpointCombination.')
    call err()
  endif
  this%combinations_ = combinations(sort(combinations,compare_combinations))
contains
  function compare_combinations(this,that) result(output)
    class(*), intent(in) :: this
    class(*), intent(in) :: that
    logical              :: output
    
    select type(this); type is(QpointCombination)
      select type(that); type is(QpointCombination)
        output = this<that
      end select
    end select
  end function
end procedure

module procedure combinations_QpointStar
  output = this%combinations_
end procedure

module procedure total_power_QpointStar
  output = this%combinations_(1)%total_power()
end procedure

module procedure equality_QpointStar_QpointStar
  if (size(this%combinations_)==size(that%combinations_)) then
    output = all(this%combinations_==that%combinations_)
  else
    output = .false.
  endif
end procedure

module procedure non_equality_QpointStar_QpointStar
  output = .not. this==that
end procedure

module procedure read_QpointStar
  type(QpointCombination), allocatable :: combinations(:)
  
  select type(this); type is(QpointStar)
    combinations = QpointCombination(input)
    this = QpointStar(combinations)
  class default
    call err()
  end select
end procedure

module procedure write_QpointStar
  select type(this); type is(QpointStar)
    output = str(this%combinations_)
  class default
    call err()
  end select
end procedure

module procedure new_QpointStar_Strings
  call this%read(input)
end procedure

module procedure new_QpointStar_StringArray
  call this%read(str(input))
end procedure

module procedure generate_stars
  integer              :: no_stars
  integer, allocatable :: star_sizes(:)
  
  integer, allocatable :: star_id(:)
  integer, allocatable :: id_in_star(:)
  
  type(QpointCombination) :: combination
  
  integer :: i,j,k,ialloc
  
  ! The trivial cases are likely to be common,
  !    so special treatment is warranted.
  if (size(combinations)==0) then
    output = [QpointStar::]
    return
  elseif (size(combinations)==1) then
    output = [QpointStar(combinations)]
    return
  endif
  
  ! Process the combinations to identify the q-point stars.
  !  - There are `no_stars` stars.
  !  - The star with id `i` contains `star_sizes(i)` combinations.
  !  - Combination `combinations(i)` is the `id_in_star(i)` combination
  !       in the star with id `star_id(i)`.
  allocate( star_sizes(size(combinations)), &
          & star_id(size(combinations)),    &
          & id_in_star(size(combinations)), &
          & stat=ialloc); call err(ialloc)
  star_id = 0
  
  ! Process the first combination.
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
    
    ! If `combinations(i)` is not in an existing star. Create a new star.
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
    allocate( output(i)%combinations_(star_sizes(i)), &
            & stat=ialloc); call err(ialloc)
  enddo
  do i=1,size(combinations)
    output(star_id(i))%combinations_(id_in_star(i)) = combinations(i)
  enddo
end procedure
end submodule

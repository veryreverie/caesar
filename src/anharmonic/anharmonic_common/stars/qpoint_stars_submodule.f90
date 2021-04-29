!> Provides the implementation of the [[QpointStars(type)]] methods.
submodule (caesar_qpoint_star_module) caesar_qpoint_stars_submodule
  use caesar_stars_module
contains

module procedure new_QpointStars
  this%power = power
  this%stars = stars
end procedure

module procedure read_QpointStars
  integer                       :: power
  type(QpointStar), allocatable :: stars(:)
  
  select type(this); type is(QpointStars)
    power = int(token(input(1), 6))
    stars = QpointStar(split_into_sections(input(2:)))
    this = QpointStars(power, stars)
  class default
    call err()
  end select
end procedure

module procedure write_QpointStars
  select type(this); type is(QpointStars)
    output = [ 'q-point stars with power = '//this%power//' :', &
             & str(this%stars, separating_line='')              ]
  class default
    call err()
  end select
end procedure

module procedure new_QpointStars_Strings
  call this%read(input)
end procedure

module procedure new_QpointStars_StringArray
  call this%read(str(input))
end procedure

module procedure generate_qpoint_stars
  type(QpointCombination), allocatable :: combinations(:)
  
  integer              :: no_stars
  integer, allocatable :: star_sizes(:)
  
  integer, allocatable :: star_id(:)
  integer, allocatable :: id_in_star(:)
  
  type(QpointCombination) :: combination
  
  integer :: i,j,k,ialloc
  
  ! Generate the q-point combinations.
  combinations = generate_qpoint_combinations(qpoints,power,conserve_momentum)
  
  if (size(combinations)==0) then
    output = [QpointStar::]
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

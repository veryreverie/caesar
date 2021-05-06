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
  type(QpointCombinations), allocatable :: combinations(:)
  
  integer :: i,ialloc
  
  combinations = generate_qpoint_combinations( qpoints,          &
                                             & max_power,        &
                                             & conserve_momentum )
  allocate(output(size(combinations)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%power = combinations(i)%power
    output(i)%stars = combinations_to_stars( combinations(i)%combinations, &
                                           & qpoint_groups                 )
  enddo
end procedure
end submodule

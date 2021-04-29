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
end submodule

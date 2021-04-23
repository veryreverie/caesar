submodule (caesar_qpoint_star_module) caesar_qpoint_star_submodule
  use caesar_anharmonic_common_module
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
    
    type(QpointPower), allocatable :: lhs(:)
    type(QpointPower), allocatable :: rhs(:)
    
    integer :: i
    
    select type(this); type is(QpointCombination)
      select type(that); type is(QpointCombination)
        lhs = this%qpoints()
        rhs = that%qpoints()
        
        if (size(lhs)/=size(rhs)) then
          call print_line(ERROR//': A QpointStar may not contain &
             &QpointCombinations containing different numbers of q-points.')
          call err()
        endif
        
        do i=1,size(lhs)
          if (lhs(i)%id()<rhs(i)%id()) then
            output = .true.
            return
          elseif (lhs(i)%id()>rhs(i)%id()) then
            output = .false.
            return
          elseif (lhs(i)%power()<rhs(i)%power()) then
            output = .true.
            return
          elseif (lhs(i)%power()>rhs(i)%power()) then
            output = .false.
            return
          elseif (lhs(i)%paired_power()<rhs(i)%paired_power()) then
            output = .true.
            return
          elseif (lhs(i)%paired_power()>rhs(i)%paired_power()) then
            output = .false.
            return
          endif
        enddo
        
        call print_line(ERROR//': A QpointStar may not contain duplicate &
           &QpointCombinations.')
        call err()
      end select
    end select
  end function
end procedure

module procedure combinations_QpointStar
  output = this%combinations_
end procedure

module procedure equality_QpointStar_QpointStar
  output = this%combinations_(1)==that%combinations_(1)
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

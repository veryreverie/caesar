submodule (caesar_qpoint_combination_module) &
   & caesar_qpoint_combination_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_QpointCombination
  this%qpoints_ = qpoints(sort(qpoints%id()))
end procedure

module procedure qpoints_QpointCombination
  output = this%qpoints_
end procedure

module procedure equality_QpointCombination_QpointCombination
  if (size(this%qpoints_)==size(that%qpoints_)) then
    output = all(this%qpoints_==that%qpoints_)
  else
    output = .false.
  endif
end procedure

module procedure non_equality_QpointCombination_QpointCombination
  output = .not. this==that
end procedure

module procedure read_QpointCombination
  type(String),      allocatable :: line(:)
  type(QpointPower), allocatable :: qpoints(:)
  
  integer :: i,ialloc
  
  select type(this); type is(QpointCombination)
    ! Splitting the input by '*' separates the q-points,
    !    but also splits q-point pairs in two.
    line = split_line(input,delimiter='*')
    
    allocate(qpoints(0), stat=ialloc); call err(ialloc)
    i = 1
    do while (i<=size(line))
      ! Check if line(i) ends in a bracket.
      if (slice(line(i),len(line(i)),len(line(i)))==')') then
        ! line(i) is a q-point on its own.
        qpoints = [qpoints, QpointPower(line(i))]
        i = i+1
      else
        ! line(i) and line(i+1) together make a q-point pair.
        qpoints = [qpoints, QpointPower(line(i)//'*'//line(i+1))]
        i = i+2
      endif
    enddo
    
    this = QpointCombination(qpoints)
  class default
    call err()
  end select
end procedure

module procedure write_QpointCombination
  select type(this); type is(QpointCombination)
    output = join(this%qpoints_, delimiter='*')
  class default
    call err()
  end select
end procedure

module procedure new_QpointCombination_String
  call this%read(input)
end procedure
end submodule

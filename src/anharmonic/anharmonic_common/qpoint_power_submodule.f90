submodule (caesar_qpoint_power_module) caesar_qpoint_power_submodule
  use caesar_anharmonic_common_module
contains

module procedure new_QpointPower
  if (present(paired_id).neqv.present(paired_power)) then
    call print_line(CODE_ERROR//': If one of paired_id and paired_power is &
      &present then the other must also be present.')
    call err()
  endif
  
  if (present(paired_id)) then
    if (paired_id==id .and. paired_power/=power) then
      call print_line(ERROR//': A QpointPower with matching id and paired_id &
         &must also have matching power and paired_power.')
      call err()
    endif
    
    if (id<=paired_id) then
      this%id_ = id
      this%power_ = power
      this%paired_id_ = paired_id
      this%paired_power_ = paired_power
    else
      this%id_ = paired_id
      this%power_ = paired_power
      this%paired_id_ = id
      this%paired_power_ = power
    endif
  else
    this%id_ = id
    this%power_ = power
    this%paired_id_ = id
    this%paired_power_ = power
  endif
end procedure

module procedure id_QpointPower
  output = this%id_
end procedure

module procedure power_QpointPower
  output = this%power_
end procedure

module procedure paired_id_QpointPower
  output = this%paired_id_
end procedure

module procedure paired_power_QpointPower
  output = this%paired_power_
end procedure

module procedure equality_QpointPower_QpointPower
  output = this%id_==that%id_               &
   & .and. this%power_==that%power_         &
   & .and. this%paired_id_==that%paired_id_ &
   & .and. this%paired_power_==that%paired_power_
end procedure

module procedure non_equality_QpointPower_QpointPower
  output = .not. this==that
end procedure

module procedure read_QpointPower
  type(String), allocatable :: line(:)
  
  integer              :: id
  integer              :: power
  integer, allocatable :: paired_id
  integer, allocatable :: paired_power
  
  select type(this); type is(QpointPower)
    line = tokens(input, delimiters=['(','q','^','*',')'])
    if (size(line)==2) then
      id = int(line(1))
      power = int(line(2))
    elseif (size(line)==4) then
      id = int(line(1))
      power = int(line(2))
      paired_id = int(line(3))
      paired_power = int(line(4))
    else
      call print_line(CODE_ERROR//': Unable to parse string to QpointPower:')
      call print_line(input)
      call err()
    endif
    this = QpointPower(id,power,paired_id,paired_power)
  class default
    call err()
  end select
end procedure

module procedure write_QpointPower
  select type(this); type is(QpointPower)
    if (this%id_==this%paired_id_) then
      output = '(q'//this%id_//'^'//this%power_//')'
    else
      output = '(q'//this%id_//'^'//this%power_//'*q'// &
             & this%paired_id_//'^'//this%paired_power_//')'
    endif
  class default
    call err()
  end select
end procedure

module procedure new_QpointPower_String
  call this%read(input)
end procedure
end submodule

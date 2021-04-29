submodule (caesar_qpoint_power_module) caesar_qpoint_power_submodule
  use caesar_stars_module
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
  
  if (this%power_<0 .or. this%paired_power_<0 .or. this%total_power()<1) then
    call print_line(CODE_ERROR//': power and paired_power must be &
       &non-negative, and total_power must be >0.')
    call err()
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

module procedure total_power_QpointPower
  if (this%paired_id_==this%id_) then
    output = this%power_
  else
    output = this%power_ + this%paired_power_
  endif
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

module procedure lt_QpointPower_QpointPower
  if (this%id_<that%id_) then
    output = .true.
  elseif (this%id_>that%id_) then
    output = .false.
  elseif (this%power_+this%paired_power_<that%power_+that%paired_power_) then
    output = .false.
  elseif (this%power_+this%paired_power_>that%power_+that%paired_power_) then
    output = .true.
  elseif (this%power_<that%power_) then
    output = .false.
  elseif (this%power_>that%power_) then
    output = .true.
  else
    output = .false.
  endif
end procedure

module procedure le_QpointPower_QpointPower
  output = this==that .or. this<that
end procedure

module procedure gt_QpointPower_QpointPower
  output = .not. this<=that
end procedure

module procedure ge_QpointPower_QpointPower
  output = .not. this<that
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

module procedure conjg_QpointPower
  output%id_ = this%id_
  output%power_ = this%paired_power_
  output%paired_id_ = this%paired_id_
  output%paired_power_ = this%power_
end procedure

module procedure operate_Group_QpointPower
  output = QpointPower( qpoint_group*this%id_,        &
                      & this%power_,                  &
                      & qpoint_group*this%paired_id_, &
                      & this%paired_power_            )
end procedure

module procedure generate_qpoint_powers
  integer :: total_power
  integer :: power
  
  if (qpoint%is_paired_qpoint()) then
    output = [( QpointPower( qpoint%id,     &
              &              total_power ), &
              & total_power=1,              &
              & max_power                   )]
  else
    output = [( ( QpointPower( qpoint%id,                     &
              &                power,                         &
              &                qpoint%paired_qpoint_id,       &
              &                total_power-power        ),    &
              &   power=0,                                    &
              &   total_power                              ), &
              & total_power=1,                                &
              & max_power                                     )]
  endif
end procedure
end submodule

submodule (caesar_group_module) caesar_group_submodule
  use caesar_algebra_module
contains

module procedure new_Group
  this%operation = operation
end procedure

module procedure size_Group
  output = size(this%operation)
end procedure

module procedure inverse_Group
  integer, allocatable :: operation(:)
  
  integer :: i,ialloc
  
  allocate(operation(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    operation(this%operation(i)) = i
  enddo
  output = Group(operation)
end procedure

module procedure operate_Group_integer
  output = this%operation(operand)
end procedure

module procedure operate_Group_Group
  if (size(this)/=size(operand)) then
    call print_line('Error: groups can only operate on other groups of the &
       & same size.')
    call err()
  endif
  
  output = Group(this*operand%operation)
end procedure

module procedure equality_Group_Group
  output = all(this%operation==that%operation)
end procedure

module procedure non_equality_Group_Group
  output = .not. this==that
end procedure

module procedure make_identity_group
  integer :: i
  
  output = Group([(i,i=1,group_size)])
end procedure

module procedure read_Group
  select type(this); type is(Group)
    this = Group(int(split_line(input)))
  end select
end procedure

module procedure write_Group
  select type(this); type is(Group)
    output = join(this%operation)
  end select
end procedure

module procedure new_Group_String
  call this%read(input)
end procedure
end submodule

module group_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  type Group
    integer, allocatable :: operation(:)
  end type
  
  interface new
    module procedure new_Group
  end interface
  
  interface drop
    module procedure drop_Group
  end interface
  
  interface assignment(=)
    module procedure assign_Group
  end interface
  
  interface size
    module procedure size_Group
  end interface
contains

subroutine new_Group(this,no_elements)
  implicit none
  
  type(Group), intent(out) :: this
  integer,     intent(in)  :: no_elements
  
  allocate(this%operation(no_elements))
end subroutine

subroutine drop_Group(this)
  implicit none
  
  type(Group), intent(inout) :: this
  
  deallocate(this%operation)
end subroutine

subroutine assign_Group(output,input)
  implicit none
  
  integer,     intent(in)  :: input(:)
  type(Group), intent(out) :: output
  
  call new(output,size(input))
  output%operation = input
end subroutine

function size_Group(this) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer                 :: output
  
  output = size(this%operation)
end function

function operate(this,input) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer,     intent(in) :: input
  integer                 :: output
  
  output = this%operation(input)
end function

function read_group_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(Group), allocatable :: this(:)
  
  type(String), allocatable :: contents(:)
  type(String), allocatable :: line(:)
  
  integer :: i
  
  contents = read_lines(filename)
  
  allocate(this(size(contents)))
  do i=1,size(contents)
    line = split(contents(i))
    call new(this(i),size(line))
    this(i)%operation = int(line)
  enddo
end function

subroutine write_group_file(this,filename)
  implicit none
  
  type(Group),  intent(in) :: this(:)
  type(String), intent(in) :: filename
  
  integer      :: group_file
  
  integer :: i
  
  group_file = open_write_file(filename)
  do i=1,size(this)
    call print_line(group_file,this(i)%operation)
  enddo
  close(group_file)
end subroutine
end module

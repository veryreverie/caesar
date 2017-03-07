module group_module
  implicit none
  
  type Group
    integer, allocatable, private :: operations(:,:)
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

subroutine new_Group(this,no_operations,no_elements)
  implicit none
  
  type(Group), intent(out) :: this
  integer,     intent(in) :: no_operations
  integer,     intent(in) :: no_elements
  
  allocate(this%operations(no_elements,no_operations))
end subroutine

subroutine drop_Group(this)
  implicit none
  
  type(Group), intent(inout) :: this
  
  deallocate(this%operations)
end subroutine

subroutine assign_Group(output,input)
  implicit none
  
  integer,     intent(in)  :: input(:,:)
  type(Group), intent(out) :: output
  
  call new(output,size(input,2),size(input,1))
  output%operations = input
end subroutine

function size_Group(this) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer                 :: output
  
  output = size(this%operations,2)
end function

function operate(input,this,operation) result(output)
  implicit none
  
  integer,     intent(in) :: input
  type(Group), intent(in) :: this
  integer,     intent(in) :: operation
  integer                 :: output
  
  output = this%operations(input,operation)
end function

function read_group_file(filename) result(this)
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: filename
  type(Group)              :: this
  
  type(String), allocatable :: contents(:)
  
  integer :: i
  
  contents = read_lines(filename)
  
  call new(this,size(contents),size(split(contents(1))))
  do i=1,size(contents)
    this%operations(:,i) = int(split(contents(i)))
  enddo
end function

subroutine write_group_file(this,filename)
  use string_module
  use file_module
  implicit none
  
  type(Group),  intent(in)  :: this
  type(String), intent(in) :: filename
  
  integer      :: group_file
  type(String) :: line
  
  integer :: i,j
  
  group_file = open_write_file(filename)
  do i=1,size(this%operations,2)
    line = ''
    do j=1,size(this%operations,1)
      line = line//' '//this%operations(j,i)
    enddo
    call print_line(group_file,line)
  enddo
  close(group_file)
end subroutine
end module

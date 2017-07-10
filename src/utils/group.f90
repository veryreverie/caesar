! ======================================================================
! A group operating on integers.
! ======================================================================

! Represents e.g. the action of a symmetry on the system.
! If a symmetry operates on a three-atom system, and takes:
!    atom 1 -> atom 2
!    atom 2 -> atom 3
!    atom 3 -> atom 1
! Then:
!    symmetry_group*1 = 2
!    symmetry_group*2 = 3
!    symmetry_group*3 = 3

module group_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! The group class.
  type Group
    integer, allocatable :: operation(:)
  end type
  
  ! Allocates the group.
  interface new
    module procedure new_Group
  end interface
  
  ! Create the group from an array of integers.
  interface assignment(=)
    module procedure assign_Group
  end interface
  
  ! Inquiries about the group.
  interface size
    module procedure size_Group
  end interface
  
  interface operator(==)
    module procedure equality_Group_Group
  end interface
  
  interface operator(/=)
    module procedure non_equality_Group_Group
  end interface
  
  ! Operate with the group on either an integer or another group.
  interface operator(*)
    module procedure operate_Group_integer
    module procedure operate_Group_Group
  end interface
contains

subroutine new_Group(this,no_elements)
  implicit none
  
  type(Group), intent(out) :: this
  integer,     intent(in)  :: no_elements
  
  allocate(this%operation(no_elements))
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

function equality_Group_Group(a,b) result(output)
  implicit none
  
  type(Group), intent(in) :: a
  type(Group), intent(in) :: b
  logical                 :: output
  
  output = all(a%operation==b%operation)
end function

function non_equality_Group_Group(a,b) result(output)
  implicit none
  
  type(Group), intent(in) :: a
  type(Group), intent(in) :: b
  logical                 :: output
  
  output = .not. a==b
end function
  
! ----------------------------------------------------------------------
! Group operations.
! ----------------------------------------------------------------------
function operate_Group_integer(this,input) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer,     intent(in) :: input
  integer                 :: output
  
  output = this%operation(input)
end function

function operate_Group_Group(a,b) result(output)
  implicit none
  
  type(Group), intent(in) :: a
  type(Group), intent(in) :: b
  type(Group)             :: output
  
  integer              :: i
  
  if (size(a)/=size(b)) then
    call print_line('Error: groups can only operate on other groups of the &
       & same size.')
    call err()
  endif
  
  call new(output, size(a))
  do i=1,size(a)
    output%operation(i) = a%operation(b%operation(i))
  enddo
end function

! ----------------------------------------------------------------------
! I/O operations with the group.
! ----------------------------------------------------------------------
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

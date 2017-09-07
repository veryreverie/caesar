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
!    symmetry_group*3 = 1
! Also:
!    (symmetry_group*symmetry_group)*1 = symmetry_group*2 = 3
module group_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: Group
  public :: size
  public :: write_group_file
  public :: read_group_file
  
  ! The group class.
  type :: Group
    integer, allocatable :: operation(:)
  contains
    generic, public  :: assignment (= ) => assign_Group
    generic, public  :: operator   (==) => equality_Group_Group
    generic, public  :: operator   (/=) => non_equality_Group_Group
    generic, public  :: operator   (* ) => operate_Group_integer, &
                                         & operate_Group_Group
    
    procedure, private :: assign_Group
    procedure, private :: equality_Group_Group
    procedure, private :: non_equality_Group_Group
    procedure, private :: operate_Group_integer
    procedure, private :: operate_Group_Group
  end type
  
  interface size
    module procedure size_Group
  end interface
  
contains

! ----------------------------------------------------------------------
! Inquiry, assignment and comparisons.
! ----------------------------------------------------------------------
! Create the group from an array of integers.
subroutine assign_Group(this,that)
  implicit none
  
  class(Group), intent(out) :: this
  integer,      intent(in)  :: that(:)
  
  this%operation = that
end subroutine

! Equality with another group.
function equality_Group_Group(this,that) result(output)
  implicit none
  
  Class(Group), intent(in) :: this
  type(Group),  intent(in) :: that
  logical                  :: output
  
  output = all(this%operation==that%operation)
end function

! Non-equality with another group.
function non_equality_Group_Group(this,that) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  type(Group),  intent(in) :: that
  logical                  :: output
  
  output = .not. this==that
end function
  
! The number of operators in the group.
function size_Group(this) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer                 :: output
  
  output = size(this%operation)
end function

! ----------------------------------------------------------------------
! Group operations.
! ----------------------------------------------------------------------
function operate_Group_integer(this,operand) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  integer,      intent(in) :: operand
  integer                  :: output
  
  output = this%operation(operand)
end function

! Defined s.t. (Group*Group)*i == Group*(Group*i) for all i.
function operate_Group_Group(this,operand) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  type(Group),  intent(in) :: operand
  type(Group)              :: output
  
  integer :: i,ialloc
  
  if (size(this)/=size(operand)) then
    call print_line('Error: groups can only operate on other groups of the &
       & same size.')
    call err()
  endif
  
  allocate(output%operation(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output%operation(i) = this%operation(operand*i)
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
  
  integer :: i
  
  contents = read_lines(filename)
  
  allocate(this(size(contents)))
  do i=1,size(contents)
    this(i)%operation = int(split(contents(i)))
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
    call print_line(group_file, this(i)%operation)
  enddo
  close(group_file)
end subroutine
end module

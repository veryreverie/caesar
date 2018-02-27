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
  
  use stringable_module
  implicit none
  
  private
  
  public :: size
  public :: make_identity_group
  
  ! The group class.
  type, public, extends(Stringable) :: Group
    integer, allocatable :: operation(:)
  contains
    generic, public  :: operator   (==) => equality_Group_Group
    generic, public  :: operator   (/=) => non_equality_Group_Group
    generic, public  :: operator   (* ) => operate_Group_integer, &
                                         & operate_Group_Group
    
    procedure, private :: equality_Group_Group
    procedure, private :: non_equality_Group_Group
    procedure, private :: operate_Group_integer
    procedure, private :: operate_Group_Group
    
    ! I/O.
    procedure :: str => str_Group
  end type
  
  interface size
    module procedure size_Group
  end interface
  
contains

! ----------------------------------------------------------------------
! Comparisons with other groups.
! ----------------------------------------------------------------------
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
  
! ----------------------------------------------------------------------
! Group properties.
! ----------------------------------------------------------------------
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
! Generates the identity group.
! ----------------------------------------------------------------------
function make_identity_group(group_size) result(output)
  implicit none
  
  integer, intent(in) :: group_size
  type(Group)         :: output
  
  integer :: i,ialloc
  
  allocate(output%operation(group_size), stat=ialloc); call err(ialloc)
  do i=1,group_size
    output%operation(i) = i
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
function str_Group(this) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  type(String)             :: output
  
  output = join(this%operation)
end function
end module

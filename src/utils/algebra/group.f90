! ======================================================================
! A group operating on integers.
! ======================================================================
! Represents e.g. the action of a symmetry on the system.
! If a symmetry operates on a three-atom system, and takes:
!    atom 1 -> atom 2
!    atom 2 -> atom 3
!    atom 3 -> atom 1
! Then:
!    group*1 = 2
!    group*2 = 3
!    group*3 = 1
! Also:
!    (group*group)*1 = group*(group*1) = group*2 = 3
module caesar_operator_group_module
  use caesar_foundations_module
  use caesar_io_module
  implicit none
  
  private
  
  public :: Group
  public :: size
  public :: operator(*)
  public :: operator(==)
  public :: operator(/=)
  public :: make_identity_group
  
  ! The group class.
  type, extends(Stringable) :: Group
    integer, allocatable :: operation(:)
  contains
    procedure, public :: inverse => inverse_Group
    ! I/O.
    procedure, public :: read  => read_Group
    procedure, public :: write => write_Group
  end type
  
  interface Group
    module procedure new_Group
    module procedure new_Group_String
  end interface
  
  interface size
    module procedure size_Group
  end interface
  
  interface operator(*)
    module procedure operate_Group_integer
    module procedure operate_Group_Group
  end interface
  
  interface operator(==)
    module procedure equality_Group_Group
  end interface
  
  interface operator(/=)
    module procedure non_equality_Group_Group
  end interface
contains
  
! ----------------------------------------------------------------------
! Constructor and size() function.
! ----------------------------------------------------------------------
function new_Group(operation) result(this)
  implicit none
  
  integer, intent(in) :: operation(:)
  type(Group)         :: this
  
  this%operation = operation
end function

function size_Group(this) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer                 :: output
  
  output = size(this%operation)
end function

! ----------------------------------------------------------------------
! Get the inverse group, i.e. the group g'(g) such that g*g' is the identity.
! ----------------------------------------------------------------------
impure elemental function inverse_Group(this) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  type(Group)              :: output
  
  integer, allocatable :: operation(:)
  
  integer :: i,ialloc
  
  allocate(operation(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    operation(this%operation(i)) = i
  enddo
  output = Group(operation)
end function

! ----------------------------------------------------------------------
! Group operations.
! ----------------------------------------------------------------------
impure elemental function operate_Group_integer(this,operand) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  integer,     intent(in) :: operand
  integer                 :: output
  
  output = this%operation(operand)
end function

! Defined s.t. (Group*Group)*i == Group*(Group*i) for all i.
impure elemental function operate_Group_Group(this,operand) result(output)
  implicit none
  
  type(Group), intent(in) :: this
  type(Group), intent(in) :: operand
  type(Group)             :: output
  
  if (size(this)/=size(operand)) then
    call print_line('Error: groups can only operate on other groups of the &
       & same size.')
    call err()
  endif
  
  output = Group(this*operand%operation)
end function

! ----------------------------------------------------------------------
! Comparisons with other groups.
! ----------------------------------------------------------------------
! Equality with another group.
impure elemental function equality_Group_Group(this,that) result(output)
  implicit none
  
  Class(Group), intent(in) :: this
  type(Group),  intent(in) :: that
  logical                  :: output
  
  output = all(this%operation==that%operation)
end function

! Non-equality with another group.
impure elemental function non_equality_Group_Group(this,that) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  type(Group),  intent(in) :: that
  logical                  :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Generates the identity group.
! ----------------------------------------------------------------------
function make_identity_group(group_size) result(output)
  implicit none
  
  integer, intent(in) :: group_size
  type(Group)         :: output
  
  integer :: i
  
  output = Group([(i,i=1,group_size)])
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_Group(this,input)
  implicit none
  
  class(Group), intent(out) :: this
  type(String), intent(in)  :: input
  
  select type(this); type is(Group)
    this = Group(int(split_line(input)))
  end select
end subroutine

function write_Group(this) result(output)
  implicit none
  
  class(Group), intent(in) :: this
  type(String)             :: output
  
  select type(this); type is(Group)
    output = join(this%operation)
  end select
end function

impure elemental function new_Group_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(Group)              :: this
  
  call this%read(input)
end function
end module

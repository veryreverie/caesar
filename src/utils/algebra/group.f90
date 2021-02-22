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
module caesar_group_module
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
    ! ----------------------------------------------------------------------
    ! Constructor and size() module function.
    ! ----------------------------------------------------------------------
    module function new_Group(operation) result(this) 
      integer, intent(in) :: operation(:)
      type(Group)         :: this
    end function
  end interface
  
  interface size
    module function size_Group(this) result(output) 
      type(Group), intent(in) :: this
      integer                 :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Get the inverse group, i.e. the group g'(g) such that g*g' is the identity.
    ! ----------------------------------------------------------------------
    impure elemental module function inverse_Group(this) result(output) 
      class(Group), intent(in) :: this
      type(Group)              :: output
    end function
  end interface
  
  interface operator(*)
    ! ----------------------------------------------------------------------
    ! Group operations.
    ! ----------------------------------------------------------------------
    impure elemental module function operate_Group_integer(this,operand) &
       & result(output) 
      type(Group), intent(in) :: this
      integer,     intent(in) :: operand
      integer                 :: output
    end function
  
    ! Defined s.t. (Group*Group)*i == Group*(Group*i) for all i.
    impure elemental module function operate_Group_Group(this,operand) &
       & result(output) 
      type(Group), intent(in) :: this
      type(Group), intent(in) :: operand
      type(Group)             :: output
    end function
  end interface
  
  interface operator(==)
    ! ----------------------------------------------------------------------
    ! Comparisons with other groups.
    ! ----------------------------------------------------------------------
    ! Equality with another group.
    impure elemental module function equality_Group_Group(this,that) &
       & result(output) 
      Class(Group), intent(in) :: this
      type(Group),  intent(in) :: that
      logical                  :: output
    end function
  end interface
  
  interface operator(/=)
    ! Non-equality with another group.
    impure elemental module function non_equality_Group_Group(this,that) &
       & result(output) 
      class(Group), intent(in) :: this
      type(Group),  intent(in) :: that
      logical                  :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates the identity group.
    ! ----------------------------------------------------------------------
    module function make_identity_group(group_size) result(output) 
      integer, intent(in) :: group_size
      type(Group)         :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_Group(this,input) 
      class(Group), intent(out) :: this
      type(String), intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_Group(this) result(output) 
      class(Group), intent(in) :: this
      type(String)             :: output
    end function
  end interface
  
  interface Group
    impure elemental module function new_Group_String(input) result(this) 
      type(String), intent(in) :: input
      type(Group)              :: this
    end function
  end interface
end module

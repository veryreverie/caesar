! ======================================================================
! A product of harmonic states along each mode in a degenerate subspace.
! ======================================================================
module subspace_product_state_module
  use common_module
  
  use single_mode_state_module
  implicit none
  
  private
  
  public :: SubspaceProductState
  public :: size
  
  type, extends(Stringable) :: SubspaceProductState
    integer, allocatable :: state_ids(:)
  contains
    procedure, public :: read  => read_SubspaceProductState
    procedure, public :: write => write_SubspaceProductState
  end type
  
  interface SubspaceProductState
    module procedure new_SubspaceProductState
    module procedure new_SubspaceProductState_String
  end interface
  
  interface size
    module procedure size_SubspaceProductState
  end interface
contains
! ----------------------------------------------------------------------
! Basic functionality:
!    - constructor.
!    - size() function.
! ----------------------------------------------------------------------
function new_SubspaceProductState(state_ids) result(this)
  implicit none
  
  integer, intent(in)        :: state_ids(:)
  type(SubspaceProductState) :: this
  
  this%state_ids = state_ids
end function

function size_SubspaceProductState(this) result(output)
  implicit none
  
  type(SubspaceProductState), intent(in) :: this
  integer                                :: output
  
  output = size(this%state_ids)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SubspaceProductState(this,input)
  implicit none
  
  class(SubspaceProductState), intent(out) :: this
  type(String),                intent(in)  :: input
  
  select type(this); type is(SubspaceProductState)
    this = SubspaceProductState(int(split_line(input)))
  end select
end subroutine

function write_SubspaceProductState(this) result(output)
  implicit none
  
  class(SubspaceProductState), intent(in) :: this
  type(String)                            :: output
  
  select type(this); type is(SubspaceProductState)
    output = join(this%state_ids)
  end select
end function

impure elemental function new_SubspaceProductState_String(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input
  type(SubspaceProductState) :: this
  
  this = input
end function
end module

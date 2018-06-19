! ======================================================================
! Provides a shared-pointer like counter.
! ======================================================================
module shared_counter_submodule
  use io_basic_module
  
  use shared_counter_bugfix_submodule
  implicit none
  
  private
  
  public :: SharedCounter
  public :: assignment(=)
  
  type :: SharedCounter
    integer, pointer, private :: no_pointers_ => null()
  contains
    final :: final_SharedCounter
    
    procedure :: is_only_pointer
    procedure :: print => print_SharedCounter
  end type
  
  interface assignment(=)
    module procedure assign_SharedCounter_SharedCounter
  end interface
  
  interface SharedCounter
    module procedure new_SharedCounter
  end interface
contains

function new_SharedCounter() result(this)
  implicit none
  
  type(SharedCounter) :: this
  
  integer :: ialloc
  
  allocate(this%no_pointers_, stat=ialloc); call err(ialloc)
  if (SHARED_COUNTER_BUG) then
    this%no_pointers_ = 0
  else
    this%no_pointers_ = 1
  endif
end function

subroutine assign_SharedCounter_SharedCounter(output,input)
  implicit none
  
  type(SharedCounter), intent(out) :: output
  type(SharedCounter), intent(in)  :: input
  
  output%no_pointers_ => input%no_pointers_
  output%no_pointers_ =  output%no_pointers_ + 1
end subroutine

impure elemental subroutine final_SharedCounter(this)
  implicit none
  
  type(SharedCounter), intent(inout) :: this
  
  if (associated(this%no_pointers_)) then
    this%no_pointers_ = this%no_pointers_ - 1
  endif
end subroutine

function is_only_pointer(this) result(output)
  implicit none
  
  class(SharedCounter), intent(in) :: this
  logical                          :: output
  
  if (.not. associated(this%no_pointers_)) then
    output = .false.
  else
    output = this%no_pointers_==1
  endif
end function

subroutine print_SharedCounter(this)
  implicit none
  
  class(SharedCounter), intent(in) :: this
  
  if (associated(this%no_pointers_)) then
    call print_line(this%no_pointers_)
  else
    call print_line('NULL')
  endif
end subroutine
end module

!> Provides a shared-pointer like counter.
module caesar_shared_counter_module
  use caesar_shared_counter_bugfix_module
  implicit none
  
  private
  
  public :: SharedCounter
  public :: assignment(=)
  
  !> A shared counter, which counts how many copies exist at once.
  type :: SharedCounter
    !> The number of copies of a given counter.
    integer, pointer, private :: no_copies_ => null()
  contains
    generic,   public  :: assignment(=) => &
                        & assign_SharedCounter_SharedCounter
    procedure, private :: assign_SharedCounter_SharedCounter
    
    final :: final_SharedCounter
    
    procedure, public :: is_only_copy => is_only_copy_SharedCounter
  end type
  
  interface SharedCounter
    !> Initialise a [[SharedCounter(type)]], and set `no_copies_` to one.
    module function new_SharedCounter() result(this) 
      type(SharedCounter) :: this
    end function
  end interface
  
  interface
    !> Copy a [[SharedCounter(type)]] from another [[SharedCounter(type)]].
    !> The copy shares a counter with the original, and this counter is
    !>    incremented when the copy happens.
    module subroutine assign_SharedCounter_SharedCounter(output,input) 
      class(SharedCounter), intent(out) :: output
      class(SharedCounter), intent(in)  :: input
    end subroutine
    
    !> Destruct a [[SharedCounter(type)]]. The counter is decremented, and if
    !>    this is the only copy the memory is freed up.
    impure elemental module subroutine final_SharedCounter(this) 
      type(SharedCounter), intent(inout) :: this
    end subroutine
    
    !> Returns `true` if this [[SharedCounter(type)]] has no other copies which
    !>    have not been destructed.
    module function is_only_copy_SharedCounter(this) result(output) 
      class(SharedCounter), intent(in) :: this
      logical                          :: output
    end function
  end interface
end module

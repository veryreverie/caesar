! ======================================================================
! Combines the Clock and MemoryChecker objects.
! ======================================================================
module profiler_module
  use precision_module
  
  use clock_module
  use memory_checker_module
  implicit none
  
  private
  
  public :: Profiler
  
  type :: Profiler
    type(Clock),         private :: clock_
    type(MemoryChecker), private :: memory_checker_
  contains
    procedure, public :: reset_clock => reset_clock_Profiler
    procedure, public :: time => time_Profiler
    procedure, public :: memory => memory_Profiler
  end type
  
  interface Profiler
    module procedure new_Profiler
  end interface
contains

function new_Profiler() result(this)
  implicit none
  
  type(Profiler) :: this
  
  this%clock_ = Clock()
  this%memory_checker_ = MemoryChecker()
end function

subroutine reset_clock_Profiler(this)
  implicit none
  
  class(Profiler), intent(inout) :: this
  
  call this%clock_%reset()
end subroutine

function time_Profiler(this) result(output)
  implicit none
  
  class(Profiler), intent(in) :: this
  real(dp)                    :: output
  
  output = this%clock_%time()
end function

function memory_Profiler(this) result(output)
  implicit none
  
  class(Profiler), intent(in) :: this
  integer                     :: output
  
  output = this%memory_checker_%memory()
end function
end module

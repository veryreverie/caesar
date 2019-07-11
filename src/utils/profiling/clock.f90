! ======================================================================
! Provides a clock for code timing purposes.
! ======================================================================
module clock_module
  use precision_module
  implicit none
  
  private
  
  public :: Clock
  
  type :: Clock
    real(dp), private :: start_time_
  contains
    procedure, public :: time
    procedure, public :: reset
  end type
  
  interface Clock
    module procedure new_Clock
  end interface
contains

impure elemental function new_Clock() result(this)
  implicit none
  
  type(Clock) :: this
  
  call cpu_time(this%start_time_)
end function

impure elemental function time(this) result(output)
  implicit none
  
  class(Clock), intent(in) :: this
  real(dp)                 :: output
  
  real(dp) :: current_time
  
  call cpu_time(current_time)
  output = current_time - this%start_time_
end function

impure elemental subroutine reset(this)
  implicit none
  
  class(Clock), intent(inout) :: this
  
  call cpu_time(this%start_time_)
end subroutine
end module

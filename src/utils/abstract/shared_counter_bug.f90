! ======================================================================
! Defines the variable SHARED_COUNTER_BUG to .true.
! ======================================================================
module shared_counter_bugfix_module
  implicit none
  
  private
  
  public :: SHARED_COUNTER_BUG
  
  logical, parameter :: SHARED_COUNTER_BUG = .true.
contains
end module

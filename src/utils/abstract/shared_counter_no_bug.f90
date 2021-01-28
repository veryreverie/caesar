! ======================================================================
! Defines the variable SHARED_COUNTER_BUG to .false.
! ======================================================================
module caesar_shared_counter_bugfix_module
  implicit none
  
  private
  
  public :: SHARED_COUNTER_BUG
  
  logical, parameter :: SHARED_COUNTER_BUG = .false.
contains
end module

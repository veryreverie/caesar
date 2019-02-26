! ======================================================================
! Provides gfortran implementations of compiler-specific functions.
! ======================================================================
module compiler_specific_module
  implicit none
  
  private
  
  public :: err_implementation
contains

! Aborts with a stacktrace.
subroutine err_implementation()
  implicit none
  
  call backtrace
  
  stop 1
end subroutine
end module

! ======================================================================
! Provides gfortran implementations of compiler-specific functions.
! ======================================================================
module compiler_specific_module
  implicit none
  
  private
  
  public :: quit_implementation
  public :: err_implementation
contains

! Aborts without a stacktrace.
subroutine quit_implementation()
  implicit none
  
  stop 1
end subroutine

! Aborts with a stacktrace.
subroutine err_implementation()
  implicit none
  
  call backtrace
  
  stop 1
end subroutine
end module

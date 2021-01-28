! ======================================================================
! Provides gfortran implementations of compiler-specific functions.
! ======================================================================
module caesar_compiler_specific_module
  implicit none
  
  private
  
  public :: quit_implementation
  public :: err_implementation
contains

! Aborts without a stacktrace.
subroutine quit_implementation()
  implicit none
  
  call exit(1)
end subroutine

! Aborts with a stacktrace.
subroutine err_implementation()
  implicit none
  
  call backtrace
  
  call exit(1)
end subroutine
end module

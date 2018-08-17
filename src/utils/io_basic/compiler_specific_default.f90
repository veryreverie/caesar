! ======================================================================
! Provides default implementations of compiler-specific functions.
! ======================================================================
module compiler_specific_submodule
  implicit none
  
  private
  
  public :: err_implementation
contains

! Aborts without a stacktrace.
subroutine err_implementation()
  implicit none
  
  integer, allocatable :: a(:)
  
  ! Attempt to cause a segfault.
  a(1) = 1
  
  ! If that fails, stop without a stacktrace.
  write(*,*)
  write(*,'(a)') "Stacktrace not implemented for this compiler."
  write(*,*)
  
  stop
end subroutine
end module

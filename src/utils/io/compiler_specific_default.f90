! ======================================================================
! Provides default implementations of compiler-specific functions.
! ======================================================================
module caesar_compiler_specific_module
  implicit none
  
  private
  
  public :: quit_implementation
  public :: err_implementation
contains

! Aborts without a stacktrace.
subroutine quit_implementation()
  stop 1
end subroutine

! Aborts with stacktrace.
subroutine err_implementation()
  integer, allocatable :: a(:)
  
  ! Attempt to cause a segfault.
  a(1) = 1
  
  ! If that fails, stop without a stacktrace.
  write(*,*)
  write(*,'(a)') "Stacktrace not implemented for this compiler."
  write(*,*)
  
  stop 1
end subroutine
end module

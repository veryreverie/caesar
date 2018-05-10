! ======================================================================
! Provides gfortran implementations of compiler-specific functions.
! ======================================================================
module compiler_specific_submodule
  implicit none
  
  private
  
  public :: err_implementation
contains

! Aborts with a stacktrace.
subroutine err_implementation()
  implicit none
  
  call abort
end subroutine
end module

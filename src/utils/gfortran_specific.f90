! ======================================================================
! This module provides gfortran-specific implementations of
!    compiler-specific functions.
! ======================================================================
module compiler_specific_module
  implicit none
contains

! Aborts with a stacktrace.
subroutine err_implementation()
  implicit none
  
  call abort
end subroutine
end module

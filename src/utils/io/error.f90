! ======================================================================
! Provides the err() subroutine, which aborts with a stacktrace.
! ======================================================================
! N.B. since different compilers handle stacktraces differently, this module
!    must wrap different files for each compiler. This is handled by the
!    separate compiler_specific_submodule files.
module error_submodule
  use terminal_submodule
  use compiler_specific_submodule
  implicit none
  
  private
  
  public :: ERROR
  public :: CODE_ERROR
  public :: WARNING
  public :: abort_with_stacktrace
  
  ! Coloured error strings.
  character(14), parameter :: ERROR = RED_ESC//'Error'//RESET_ESC
  character(19), parameter :: CODE_ERROR = RED_ESC//'Code Error'//RESET_ESC
  character(16), parameter :: WARNING = LIGHT_MAGENTA_ESC//'Warning'//RESET_ESC
contains

! Abort with a stacktrace.
subroutine abort_with_stacktrace()
  implicit none
  
  call err_implementation()
end subroutine
end module


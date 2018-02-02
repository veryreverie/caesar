! ======================================================================
! Provides the err() subroutine, which aborts with a stacktrace.
! ======================================================================
! N.B. since different compilers handle stacktraces differently, this module
!    must wrap different files for each compiler. This is handled by the
!    separate compiler_specific_module files.
! N.B. Everything here is re-exported in io_module, to save from having to
!    import this module explicitly.
module error_module
  use terminal_module
  implicit none
  
  ! Coloured error strings.
  character(14), parameter :: ERROR = RED_ESC//'Error'//RESET_ESC
  character(19), parameter :: CODE_ERROR = RED_ESC//'Code Error'//RESET_ESC
  character(16), parameter :: WARNING = LIGHT_MAGENTA_ESC//'Warning'//RESET_ESC
contains

! Abort with a stacktrace.
subroutine abort_with_stacktrace()
  use compiler_specific_module
  implicit none
  
  call err_implementation()
end subroutine
end module


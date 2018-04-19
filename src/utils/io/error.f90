! ======================================================================
! Provides the err() subroutine, which aborts with a stacktrace.
! Also provides coloured "Error" and "Warning" strings for messages.
! ======================================================================
module error_submodule
  use terminal_submodule
  use compiler_specific_submodule
  implicit none
  
  private
  
  public :: ERROR
  public :: CODE_ERROR
  public :: WARNING
  public :: set_error_strings_coloured
  public :: set_error_strings_uncoloured
  public :: abort_with_stacktrace
  
  ! Coloured error strings.
  character(:), allocatable, protected :: ERROR
  character(:), allocatable, protected :: CODE_ERROR
  character(:), allocatable, protected :: WARNING
contains

! Set whether or not the error strings are enclosed by
!    terminal colour escape characters.
subroutine set_error_strings_coloured()
  implicit none
  
  ERROR = RED_ESC//'Error'//RESET_ESC
  CODE_ERROR = RED_ESC//'Code Error'//RESET_ESC
  WARNING = LIGHT_MAGENTA_ESC//'Warning'//RESET_ESC
end subroutine

subroutine set_error_strings_uncoloured()
  implicit none
  
  ERROR = 'Error'
  CODE_ERROR = 'Code Error'
  WARNING = 'Warning'
end subroutine

! Abort with a stacktrace.
subroutine abort_with_stacktrace()
  implicit none
  
  call err_implementation()
end subroutine
end module


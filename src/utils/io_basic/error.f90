! ======================================================================
! Provides the err() subroutine, which aborts with a stacktrace.
! Also provides coloured "Error" and "Warning" strings for messages.
! ======================================================================
module error_module
  use terminal_module
  use compiler_specific_module
  implicit none
  
  private
  
  public :: ERROR
  public :: CODE_ERROR
  public :: WARNING
  public :: quit
  public :: err
  public :: set_error_strings_coloured
  public :: set_error_strings_uncoloured
  
  ! Coloured error strings.
  character(:), allocatable, protected :: ERROR
  character(:), allocatable, protected :: CODE_ERROR
  character(:), allocatable, protected :: WARNING
  
  interface err
    module procedure err_none
    module procedure err_allocate_flag
  end interface
contains

! ----------------------------------------------------------------------
! Aborts without a stacktrace.
! ----------------------------------------------------------------------
subroutine quit()
  implicit none
  
  call quit_implementation()
end subroutine

! ----------------------------------------------------------------------
! Aborts with a stacktrace.
! ----------------------------------------------------------------------
! Always aborts.
subroutine err_none()
  implicit none
  
  write(*,'(a)') ''
  write(*,'(a)') LIGHT_MAGENTA_ESC//'Intentionally aborting with &
                 &stacktrace.'//RESET_ESC
  write(*,'(a)') ''
  write(*,'(a)') LIGHT_MAGENTA_ESC//'compiler_specific.f90 and error.f90 &
                 &can be ignored in stacktrace.'//RESET_ESC
  write(*,'(a)') ''
  
  call err_implementation()
end subroutine

! Aborts if integer input /= 0.
! Designed for use with allocate 'stat=ierr' flags.
subroutine err_allocate_flag(this)
  implicit none
  
  integer, intent(in) :: this
  
  if (this/=0) then
    write(*,'(a)') ''
    write(*,'(a)') ERROR//': Allocation error.'
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! Set whether or not the error strings are enclosed by
!    terminal colour escape characters.
! ----------------------------------------------------------------------
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
end module


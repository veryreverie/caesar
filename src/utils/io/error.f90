!> Provides the [[err]] subroutine, which aborts with a stacktrace.
!> Also provides coloured `Error` and `Warning` strings for messages.
module caesar_error_module
  use caesar_terminal_module
  use caesar_compiler_specific_module
  implicit none
  
  private
  
  public :: ERROR
  public :: CODE_ERROR
  public :: WARNING
  public :: quit
  public :: err
  public :: set_error_strings_coloured
  public :: set_error_strings_uncoloured
  
  character(5),  target, private :: ERROR_UNCOLOURED = 'Error'
  character(10), target, private :: CODE_ERROR_UNCOLOURED = 'Code Error'
  character(7),  target, private :: WARNING_UNCOLOURED = 'Warning'
  
  character(14), target, private :: ERROR_COLOURED = &
                                  & RED_ESC//'Error'//RESET_ESC
  character(19), target, private :: CODE_ERROR_COLOURED = &
                                  & RED_ESC//'Code Error'//RESET_ESC
  character(16), target, private :: WARNING_COLOURED = &
                                  & LIGHT_MAGENTA_ESC//'Warning'//RESET_ESC
  
  character(:), pointer, protected :: ERROR => ERROR_COLOURED
  character(:), pointer, protected :: CODE_ERROR => CODE_ERROR_UNCOLOURED
  character(:), pointer, protected :: WARNING => WARNING_UNCOLOURED
  
  !!> The string `ERROR`, in red if colour enabled.
  !character(:), allocatable, protected :: ERROR
  !!> The string `CODE ERROR`, in red if colour enabled.
  !character(:), allocatable, protected :: CODE_ERROR
  !!> The string `WARNING`, in magenta if colour enabled.
  !character(:), allocatable, protected :: WARNING
  
  interface
    !> Aborts without a stacktrace.
    module subroutine quit()
    end subroutine
  end interface
    
  interface err
    !> Aborts with a stacktrace.
    module subroutine err_none()
    end subroutine

    !> Aborts with a stacktrace if `input` $\neq 0$.
    !> Designed for use with `allocate(foo, stat=ialloc)` flags.
    module subroutine err_allocate_flag(this)
      integer, intent(in) :: this
    end subroutine
  end interface

  interface
    !> Sets `ERROR`, `CODE_ERROR` and `WARNING` to be coloured.
    module subroutine set_error_strings_coloured()
    end subroutine

    !> Sets `ERROR`, `CODE_ERROR` and `WARNING` to be uncoloured.
    module subroutine set_error_strings_uncoloured()
    end subroutine
  end interface
end module

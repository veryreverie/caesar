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
  
  !> The string `Error` without terminal escape characters for colour.
  character(5),  target, private :: ERROR_UNCOLOURED = 'Error'
  !> The string `Code Error` without terminal escape characters for colour.
  character(10), target, private :: CODE_ERROR_UNCOLOURED = 'Code Error'
  !> The string `Warning` without terminal escape characters for colour.
  character(7),  target, private :: WARNING_UNCOLOURED = 'Warning'
  
  !> The string `Error` with terminal escape characters for red.
  character(14), target, private :: ERROR_COLOURED = &
                                  & RED_ESC//'Error'//RESET_ESC
  !> The string `Code Error` with terminal escape characters for red.
  character(19), target, private :: CODE_ERROR_COLOURED = &
                                  & RED_ESC//'Code Error'//RESET_ESC
  !> The string `Warning` with terminal escape characters for magenta.
  character(16), target, private :: WARNING_COLOURED = &
                                  & LIGHT_MAGENTA_ESC//'Warning'//RESET_ESC
  
  !> The string `Error`, with or without terminal escape characters for colour
  !>    as defined by [[set_error_strings_coloured]] and
  !>    [[set_error_strings_uncoloured]].
  character(:), pointer, protected :: ERROR => ERROR_COLOURED
  !> The string `Code Error`, with or without terminal escape characters for
  !>    colour as defined by [[set_error_strings_coloured]] and
  !>    [[set_error_strings_uncoloured]].
  character(:), pointer, protected :: CODE_ERROR => CODE_ERROR_COLOURED
  !> The string `Warning`, with or without terminal escape characters for
  !>    colour as defined by [[set_error_strings_coloured]] and
  !>    [[set_error_strings_uncoloured]].
  character(:), pointer, protected :: WARNING => WARNING_COLOURED
  
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

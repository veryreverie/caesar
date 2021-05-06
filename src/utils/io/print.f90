!> Provides functionallity for writing to the terminal.
module caesar_print_module
  use iso_fortran_env, only : OUTPUT_UNIT
  use caesar_terminal_module
  use caesar_error_module
  use caesar_string_base_module
  use caesar_string_module
  use caesar_print_settings_module
  implicit none
  
  private
  
  public :: print_line
  public :: print_lines
  public :: colour
  public :: set_output_unit
  public :: unset_output_unit
  
  interface print_line
    !> Prints a character or string variable, either to the terminal or
    !>    to the output file specified by OUTPUT_FILE_UNIT.
    !> Acts as the intrinsic write(*,*), but with the ability to be
    !>    redirected to a file. Provides error checking and formatting
    !>   dependent on whether the output is the terminal or a file.
    module subroutine print_line_character(line,settings)
      character(*),        intent(in)           :: line
      !> Settings determining how the print behaves.
      !> Defaults to the global [[PrintSettings(type)]]
      !>    set by `set_print_settings`.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine

    !> Prints a character or string variable, either to the terminal or
    !>    to the output file specified by OUTPUT_FILE_UNIT.
    !> Acts as the intrinsic write(*,*), but with the ability to be
    !>    redirected to a file. Provides error checking and formatting
    !>   dependent on whether the output is the terminal or a file.
    module subroutine print_line_String(line,settings)
      type(String),        intent(in)           :: line
      !> Settings determining how the print behaves.
      !> Defaults to the global [[PrintSettings(type)]]
      !>    set by `set_print_settings`.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine
  end interface
  
  interface print_lines
    !> As print_line, but for printing over multiple lines.
    module subroutine print_lines_Strings_String(lines,separating_line, &
       & settings)
      type(String),        intent(in)           :: lines(:)
      !> If present, `separating_line` will be printed between each line.
      type(String),        intent(in), optional :: separating_line
      !> Settings determining how the print behaves.
      !> Defaults to the global [[PrintSettings(type)]]
      !>    set by `set_print_settings`.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine
    
    !> As print_line, but for printing over multiple lines.
    module subroutine print_lines_Strings_character(lines,separating_line, &
       & settings)
      type(String),        intent(in)           :: lines(:)
      !> If present, `separating_line` will be printed between each line.
      character(*),        intent(in)           :: separating_line
      !> Settings determining how the print behaves.
      !> Defaults to the global [[PrintSettings(type)]]
      !>    set by `set_print_settings`.
      type(PrintSettings), intent(in), optional :: settings
    end subroutine
  end interface
  
  interface colour
    !> Adds terminal escape characters so that a string is coloured when
    !>    printed to the terminal.
    !> Does nothing if OUTPUT_FILE is set.
    module function colour_character_character(input,colour_name) &
       & result(output)
      character(*), intent(in) :: input
      character(*), intent(in) :: colour_name
      type(String)             :: output
    end function

    !> Adds terminal escape characters so that a string is coloured when
    !>    printed to the terminal.
    !> Does nothing if OUTPUT_FILE is set.
    module function colour_String_character(input,colour_name) result(output)
      character(*), intent(in) :: input
      type(String), intent(in) :: colour_name
      type(String)             :: output
    end function

    !> Adds terminal escape characters so that a string is coloured when
    !>    printed to the terminal.
    !> Does nothing if OUTPUT_FILE is set.
    module function colour_character_String(input,colour_name) result(output)
      type(String), intent(in) :: input
      character(*), intent(in) :: colour_name
      type(String)             :: output
    end function

    !> Adds terminal escape characters so that a string is coloured when
    !>    printed to the terminal.
    !> Does nothing if OUTPUT_FILE is set.
    impure elemental module function colour_String_String(input,colour_name) &
       & result(output)
      type(String), intent(in) :: input
      type(String), intent(in) :: colour_name
      type(String)             :: output
    end function
  end interface
  
  interface
    !> C tput cols interface.
    function get_terminal_width_c(width) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      integer(kind=c_int), intent(out) :: width
      logical(kind=c_bool)             :: success
    end function
  end interface
  
  interface set_output_unit
    !> Set OUTPUT_FILE_UNIT to a file unit.
    module subroutine set_output_unit(file_unit)
      integer, intent(in) :: file_unit
    end subroutine
  end interface
  
  interface unset_output_unit
    !> Unset OUTPUT_FILE_UNIT.
    module subroutine unset_output_unit()
    end subroutine
  end interface
end module

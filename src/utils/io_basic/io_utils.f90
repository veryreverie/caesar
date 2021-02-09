!> Various miscellaneous I/O operations.
module caesar_io_utils_module
  use iso_fortran_env, only : INPUT_UNIT
  use caesar_precision_module
  
  use caesar_terminal_module
  use caesar_error_module
  use caesar_string_base_module
  use caesar_string_module
  use caesar_print_module
  use caesar_intrinsics_module
  implicit none
  
  private
  
  public :: startup_io_utils
  public :: CommandLineFlag
  public :: file_exists
  public :: command_line_args
  public :: mkdir
  public :: set_working_directory
  public :: system_call
  public :: get_flag
  public :: read_line_from_user
  public :: format_path
  public :: execute_old_code
  public :: execute_python
  public :: call_caesar
  public :: parse_c_string
  
  !> Stores a flag from the command line input, and its argument.
  type :: CommandLineFlag
    character(1) :: flag
    type(String) :: argument
  end type
  
  interface file_exists
    !> Checks if a file exists, given its filename.
    module function file_exists_character(filename) result(output) 
      character(*), intent(in) :: filename
      logical                  :: output
    end function

    !> Checks if a file exists, given its filename.
    module function file_exists_string(filename) result(output) 
      type(String), intent(in) :: filename
      logical                  :: output
    end function
  end interface
  
  interface format_path
    !> Takes a directory name, and converts it into an absolute path
    !>    in standard format (without a trailing `/`).
    !> If `working_directory` is specified, then relative paths are taken to be
    !>    relative to `working_directory`. Otherwise, they are taken to be
    !>    relative to Caesar's current working directory.
    module function format_path_character(path,working_directory) &
       & result(output) 
      character(*), intent(in)           :: path
      type(String), intent(in), optional :: working_directory
      type(String)                       :: output
    end function

    !> Takes a directory name, and converts it into an absolute path
    !>    in standard format (without a trailing `/`).
    !> If `working_directory` is specified, then relative paths are taken to be
    !>    relative to `working_directory`. Otherwise, they are taken to be
    !>    relative to Caesar's current working directory.
    module function format_path_String(path,working_directory) result(output) 
      type(String), intent(in)           :: path
      type(String), intent(in), optional :: working_directory
      type(String)                       :: output
    end function
  end interface
  
  interface set_working_directory
    !> Sets the current working directory to `working_directory`.
    module subroutine set_working_directory_character(working_directory) 
      character(*), intent(in) :: working_directory
    end subroutine

    !> Sets the current working directory to `working_directory`.
    module subroutine set_working_directory_String(working_directory) 
      type(String), intent(in) :: working_directory
    end subroutine
  end interface
  
  interface call_caesar
    !> Calls `caesar` on the command line, with the given `arguments`.
    module subroutine call_caesar_character(arguments) 
      character(*), intent(in), optional :: arguments
    end subroutine

    !> Calls `caesar` on the command line, with the given `arguments`.
    module subroutine call_caesar_String(arguments) 
      type(String), intent(in) :: arguments
    end subroutine
  end interface
  
  interface command_line_args
    !> Returns an array containing the command line arguments.
    module function command_line_args() result(args) 
      type(String), allocatable :: args(:)
    end function
  end interface

  interface mkdir
    !> Make a directory, if it doesn't already exist.
    module subroutine mkdir(dirname) 
      type(String), intent(in) :: dirname
    end subroutine
  end interface

  interface system_call
    !> Makes a system call using `system.c`.
    module function system_call(this) result(result_code) 
      type(String), intent(in) :: this
      integer                  :: result_code
    end function
  end interface

  interface get_flag
    !> Reads flags from the command line via system.c.
    !> For none-flag arguments (those not preceded by `-` or `--`),
    !>    `output%flag` is set to a space, `' '`.
    !> For flags which accept arguments but do not have them,
    !>    `output%argument` is set to an empty string, `''`.
    !> For flags which do not accept arguments, `output%argument` is undefined.
    !> Aborts if an unexpected flag is passed.
    module function get_flag(args,flags_without_arguments, &
       & flags_with_arguments) result(output) 
      type(String), intent(in) :: args(:)
      type(String), intent(in) :: flags_without_arguments
      type(String), intent(in) :: flags_with_arguments
      type(CommandLineFlag)    :: output
    end function
  end interface

  interface read_line_from_user
    !> Reads a line from `stdin`.
    module function read_line_from_user() result(line) 
      type(String) :: line
    end function
  end interface

  interface get_home_directory
    !> Sets `HOME` to the path to the directory '~'.
    module function get_home_directory() result(output) 
      type(String) :: output
    end function
  end interface

  interface get_current_directory
    !> Sets `CWD` to the path to the current working directory.
    module function get_current_directory() result(output) 
      type(String) :: output
    end function
  end interface

  interface get_exe_location
    !> Sets `EXE_LOCATION` to the path to the `caesar` executable.
    module function get_exe_location() result(output) 
      type(String) :: output
    end function
  end interface

  interface startup_io_utils
    !> Initialises this module's global variables.
    !> Must be called before any other procedures in this module.
    module subroutine startup_io_utils() 
    end subroutine
  end interface
  
  interface execute_old_code
    !> Executes one of the old caesar executables.
    !> Adds the path to old Caesar to PATH in a subshell.
    !> Moves to the working directory.
    !> Runs the given file.
    module subroutine execute_old_code(filename) 
      type(String), intent(in) :: filename
    end subroutine
  end interface

  interface execute_python
    !> Executes one of the Python scripts.
    module subroutine execute_python(filename,python_path,python_arguments) 
      type(String), intent(in)           :: filename
      type(String), intent(in)           :: python_path
      type(String), intent(in), optional :: python_arguments(:)
    end subroutine
  end interface

  interface parse_c_string
    !> Converts a `NUL`-terminated C string to a Fortran [[String(type)]].
    module function parse_c_string(input) result(output) 
      character(*), intent(in) :: input
      type(String)             :: output
    end function
  end interface

  interface escape_bash_characters
    !> Adds escape characters to a string so that it can be used on the
    !>    command line.
    recursive module function escape_bash_characters(input) result(output) 
      character(*), intent(in) :: input
      type(String)             :: output
    end function
  end interface
  
  !> Interface to C functions in `system.c`.
  interface
    !> `getcwd` interface.
    function get_cwd_c(result_size, cwd) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      integer(kind=c_int),    intent(in)  :: result_size
      character(kind=c_char), intent(out) :: cwd(*)
      logical(kind=c_bool)                :: success
    end function
  
    !> `getenv('HOME')` interface.
    function get_home_c(home) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(out) :: home(*)
      logical(kind=c_bool)                :: success
    end function
  
    !> `readlink('/proc/self/exe')` interface.
    function get_exe_location_c(result_size,exe_location) bind(c) &
       & result(success)
      use, intrinsic :: iso_c_binding
      integer(kind=c_int),    intent(in)  :: result_size
      character(kind=c_char), intent(out) :: exe_location
      logical(kind=c_bool)                :: success
    end function
  
    !> Flag parsing interface.
    function get_flag_c(argc,argvs,options,flag,output) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      integer(kind=c_int),    intent(in)    :: argc
      character(kind=c_char), intent(inout) :: argvs(*)
      character(kind=c_char), intent(in)    :: options(*)
      character(kind=c_char), intent(out)   :: flag
      character(kind=c_char), intent(out)   :: output(*)
      logical(kind=c_bool)                  :: success
    end function
  
    !> System call interface.
    function system_c(input) bind(c) result(return_code)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(in) :: input(*)
      integer(kind=c_int)                :: return_code
    end function
  
    !> `mkdir` interface.
    function mkdir_c(input) bind(c) result(return_code)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(in) :: input(*)
      integer(kind=c_int)                :: return_code
    end function
  end interface
end module

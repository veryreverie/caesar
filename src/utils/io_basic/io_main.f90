! ======================================================================
! I/O operations.
! ======================================================================
module io_submodule
  use iso_fortran_env, only : INPUT_UNIT
  use precision_module
  
  use terminal_submodule
  use error_submodule
  use string_submodule
  use print_submodule
  use intrinsics_submodule
  implicit none
  
  private
  
  ! Public types.
  public :: CommandLineFlag
  
  ! IO operations.
  public :: file_exists         ! checks if a file exists
  public :: command_line_args   ! Get arguments from the command line.
  public :: mkdir               ! Make a directory.
  public :: set_io_settings     ! Sets I/O settings.
  public :: system_call         ! Makes a system call.
  public :: get_flag            ! Reads a flag from the command line.
  public :: read_line_from_user ! Reads a line from the terminal.
  public :: print_line          ! write(*,'(a)')
  public :: format_path         ! Converts any path into an absolute path.
  public :: execute_old_code    ! Runs one of the old caesar codes.
  public :: execute_python      ! Runs one of the python scripts.
  
  ! I/O settings, specifying various input/output properties.
  ! Set by set_io_settings.
  type(String) :: HOME
  type(String) :: CWD
  type(String) :: OLD_PATH
  type(String) :: PYTHON_SCRIPTS_PATH
  
  ! Command line flag and argument.
  type :: CommandLineFlag
    character(1) :: flag
    type(String) :: argument
  end type
  
  interface file_exists
    module procedure file_exists_character
    module procedure file_exists_string
  end interface
  
  interface print_line
    module procedure print_line_integer
    module procedure print_line_real
    module procedure print_line_logical
    module procedure print_line_complex
    
    module procedure print_line_integers
    module procedure print_line_reals
    module procedure print_line_logicals
    module procedure print_line_complexes
  end interface
  
  interface format_path
    module procedure format_path_character
    module procedure format_path_String
  end interface
  
  ! C system call interface.
  interface
    function system_c(input) bind(c) result(return_code)
      use, intrinsic :: iso_c_binding
      implicit none
      
      character(kind=c_char), intent(in) :: input(*)
      integer(kind=c_int)                :: return_code
    end function
  end interface
  
  ! C getcwd call interface.
  interface
    function get_cwd_c(result_size, cwd) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int),    intent(in)  :: result_size
      character(kind=c_char), intent(out) :: cwd(*)
      logical(kind=c_bool)                :: success
    end function
  end interface
  
  ! C getenv('HOME') call interface.
  interface
    function get_home_c(home) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      character(kind=c_char), intent(out) :: home(*)
      logical(kind=c_bool)                :: success
    end function
  end interface
  
  ! C readlink('/proc/self/exe') call interface.
  interface
    function get_exe_location_c(result_size,exe_location) bind(c) &
       & result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int),    intent(in)  :: result_size
      character(kind=c_char), intent(out) :: exe_location
      logical(kind=c_bool)                :: success
    end function
  end interface
  
  ! C flag parsing interface.
  interface
    function get_flag_c(argc,argvs,options,flag,output) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int),    intent(in)    :: argc
      character(kind=c_char), intent(inout) :: argvs(*)
      character(kind=c_char), intent(in)    :: options(*)
      character(kind=c_char), intent(out)   :: flag
      character(kind=c_char), intent(out)   :: output(*)
      logical(kind=c_bool)                  :: success
    end function
  end interface
contains

! ----------------------------------------------------------------------
! String manipulation.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! I/O operations.
! ----------------------------------------------------------------------

! Checks if a file exists
function file_exists_character(filename) result(output)
  implicit none
  
  character(*), intent(in) :: filename
  logical                  :: output
  
  inquire(file=filename, exist=output)
end function

function file_exists_string(filename) result(output)
  implicit none
  
  type(String), intent(in) :: filename
  logical                  :: output
  
  output = file_exists(char(filename))
end function

! Returns an array containing the command line arguments.
function command_line_args() result(args)
  implicit none

  type(String), allocatable :: args(:)
  
  ! Temporary variables.
  integer         :: i,ialloc
  integer         :: arg_count
  character(1000) :: temp_char

  arg_count = command_argument_count()
  allocate (args(arg_count+1), stat=ialloc); call err(ialloc)
  do i=0,arg_count
    call get_command_argument(i, temp_char)
    args(i+1) = trim(temp_char)
  enddo
end function

! Make a directory, if it doesn't already exist.
subroutine mkdir(dirname)
  implicit none
  
  type(String), intent(in) :: dirname
  
  integer :: result_code
  
  result_code = system_call( &
     & 'if [ ! -e '//dirname//' ]; then mkdir '//dirname//'; fi')
  if (result_code/=0) then
    call print_line('Error: failed to make directory: '//dirname)
    call err()
  endif
end subroutine

! Makes a system call via system.c.
function system_call(this) result(result_code)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: result_code
  
  result_code = system_c(char(this)//char(0))
end function

! ----------------------------------------------------------------------
! Reads flags from the command line via system.c.
! ----------------------------------------------------------------------
! For none-flag arguments (those not preceded by '-' or '--'), output%flag = ' '.
! For flags which accept arguments but do not have them, output%argument = ''.
! For flags which do not accept arguments, output%argument is undefined.
! Aborts if an unexpected flag is passed.
function get_flag(args,flags_without_arguments,flags_with_arguments) &
   & result(output)
  implicit none
  
  type(String), intent(in) :: args(:)
  type(String), intent(in) :: flags_without_arguments
  type(String), intent(in) :: flags_with_arguments
  type(CommandLineFlag)    :: output
  
  ! Arguments passed to C.
  character(1000) :: argvs
  type(String)    :: options
  character(1)    :: flag
  character(1000) :: argument
  logical         :: success
  
  ! Temporary variables.
  integer :: i
  character(1000) :: flags_char
  
  ! Convert command-line arguments into C-friendly format.
  argvs = char(join(args,char(0))//char(0))
  
  ! Convert flags into C-friendly format.
  options = '-:-:'//flags_without_arguments
  flags_char = char(flags_with_arguments)
  do i=1,len(flags_with_arguments)
    options = options//flags_char(i:i)//':'
  enddo
  options = options//char(0)
  
  ! Call getopt.
  success = get_flag_c( size(args),    &
                      & argvs,         &
                      & char(options), &
                      & flag,          &
                      & argument)
  
  ! Return result.
  if (.not. success) then
    ! All arguments have been processed.
    output%flag = ' '
    output%argument = ''
  else
    output%flag = flag
    do i=1,len(argument)
      if (argument(i:i)==char(0)) then
        output%argument = argument(:i-1)
        exit
      endif
    enddo
    
    ! The flag is unexpected.
    if (output%flag=='?') then
      call print_line('Error: unexpected flag "'//output%argument//'".')
      stop
    
    ! The flag has no argument.
    elseif (output%flag==':') then
      output%flag = output%argument
      output%argument = ''
    endif
  endif
end function

function read_line_from_user() result(line)
  implicit none
  
  type(String) :: line
  
  character(1000) :: char_line
  
  read(INPUT_UNIT,'(a)') char_line
  line = trim(char_line)
end function

! ----------------------------------------------------------------------
! Sets HOME, CWD and TERMINAL_WIDTH.
! Calls the system to find these variables.
! ----------------------------------------------------------------------
function get_home_directory() result(output)
  implicit none
  
  type(String) :: output
  
  integer, parameter  :: result_size = 1024
  
  character(result_size) :: home_dir
  logical                :: success
  integer                :: i
  
  success = get_home_c(home_dir)
  
  if (.not. success) then
    call print_line('Error: getenv("HOME") failed.')
    call err()
  endif
  
  output = ''
  do i=1,result_size
    if (home_dir(i:i)==char(0)) then
      output = home_dir(:i-1)
      exit
    endif
  enddo
  
  if (output=='') then
    call print_line('Error: home directory string not nul-terminated: '// &
       & home_dir)
    call err()
  endif
end function

function get_current_directory() result(output)
  implicit none
  
  type(String) :: output
  
  integer, parameter :: result_size = 1024
  
  character(result_size) :: current_dir
  logical                :: success
  integer                :: i
  
  success = get_cwd_c(result_size, current_dir)
  
  if (.not. success) then
    call print_line('Error: getcwd failed.')
    call err()
  endif
  
  output = ''
  do i=1,result_size
    if (current_dir(i:i)==char(0)) then
      output = current_dir(:i-1)
      exit
    endif
  enddo
  
  if (output=='') then
    call print_line('Error: cwd string not nul-terminated: '//current_dir)
    call err()
  endif
end function

function get_exe_location() result(output)
  implicit none
  
  type(String) :: output
  
  integer, parameter :: result_size = 1024
  
  character(result_size) :: exe_location
  logical                :: success
  integer                :: i
  
  success = get_exe_location_c(result_size, exe_location)
  
  if (.not. success) then
    call print_line('Error: readlink("/proc/self/exe") failed.')
    call err()
  endif
  
  output = ''
  do i=1,result_size
    if (exe_location(i:i)==char(0)) then
      output = exe_location(:i-1)
      exit
    endif
  enddo
  
  if (output=='') then
    call print_line('Error: readlink string not nul-terminated: '// &
       & exe_location)
    call err()
  endif
end function

subroutine set_io_settings()
  implicit none
  
  type(String) :: exe_location
  type(String) :: caesar_dir
  
  call set_terminal_width()
  HOME = get_home_directory()
  CWD = get_current_directory()
  exe_location = get_exe_location()
  if (  slice(exe_location,len(exe_location)-10,len(exe_location)) &
   & /= '/bin/caesar') then
    call print_line('Caesar executable in unexpected location: '//exe_location)
    call err()
  endif
  caesar_dir = slice(exe_location,1,len(exe_location)-11)
  OLD_PATH = caesar_dir//'/old'
  PYTHON_SCRIPTS_PATH = caesar_dir//'/python'
  call unset_output_unit()
end subroutine

! ----------------------------------------------------------------------
! Provides the print_line function for non-character types.
! ----------------------------------------------------------------------
subroutine print_line_integer(this)
  implicit none
  
  integer, intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_real(this)
  implicit none
  
  real(dp), intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_logical(this)
  implicit none
  
  logical, intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_complex(this)
  implicit none
  
  complex(dp), intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_integers(this)
  implicit none
  
  integer, intent(in) :: this(:)
  
  call print_line(join(this))
end subroutine

subroutine print_line_reals(this)
  implicit none
  
  real(dp), intent(in) :: this(:)
  
  call print_line(join(this))
end subroutine

subroutine print_line_logicals(this)
  implicit none
  
  logical, intent(in) :: this(:)
  
  call print_line(join(this))
end subroutine

subroutine print_line_complexes(this)
  implicit none
  
  complex(dp), intent(in) :: this(:)
  
  call print_line(join(this))
end subroutine

! ----------------------------------------------------------------------
! Takes a directory name, and converts it into an absolute path
!    in standard format (without a trailing '/').
! ----------------------------------------------------------------------
function format_path_character(path) result(output)
  implicit none
  
  character(*), intent(in) :: path
  type(String)             :: output
  
  integer :: length
  
  length = len(path)
  
  if (length==0) then
    call print_line('Error: no path provided.')
    call err()
  endif
  
  ! foo/ -> foo
  if (length>1 .and. path(length:)=='/') then
    length = length - 1
  endif
  
  if (length==1) then
    if (path(:1)=='/') then
      ! / -> /
      output = path(:1)
    elseif (path(:1)=='.') then
      ! . -> /current/working/directory
      output = CWD
    elseif (path(:1)=='~') then
      ! ~ -> /home/directory
      output = HOME
    else
      ! foo -> /current/working/directory/foo
      output = CWD//'/'//path(:1)
    endif
  else
    if (path(:1)=='/') then
      ! /foo -> /foo
      output = path(:length)
    elseif (path(:2)=='./') then
      ! ./foo -> /current/working/directory/foo
      output = CWD//path(2:length)
    elseif (path(:2)=='~/') then
      ! ~/foo -> /home/directory/foo
      output = HOME//path(2:length)
    else
      ! foo -> /current/working/directory/foo
      output = CWD//'/'//path(:length)
    endif
  endif
end function

function format_path_String(path) result(output)
  implicit none
  
  type(String), intent(in) :: path
  type(String)             :: output
  
  output = format_path(char(path))
end function

! ----------------------------------------------------------------------
! Executes one of the old caesar executables.
! ----------------------------------------------------------------------
! Adds the path to old Caesar to PATH in a subshell.
! Moves to the working directory.
! Runs the given file.
subroutine execute_old_code(wd, filename)
  implicit none
  
  type(String), intent(in) :: wd
  type(String), intent(in) :: filename
  
  integer :: result_code
  
  result_code = system_call( &
     & '( PATH='//OLD_PATH//':$PATH; cd '//wd//'; '//filename//' )')
  
  if (result_code/=0) then
    call print_line('Error: '//filename//' failed.')
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! Executes one of the Python scripts.
! ----------------------------------------------------------------------
subroutine execute_python(wd, filename, python_path)
  implicit none
  
  type(String), intent(in) :: wd
  type(String), intent(in) :: filename
  type(String), intent(in) :: python_path
  
  integer :: result_code
  
  result_code = system_call( 'cd '//wd//'; '//  &
                           & python_path//' '// &
                           & PYTHON_SCRIPTS_PATH//'/'//filename)
  
  if (result_code/=0) then
    call print_line('Error: '//filename//' failed.')
    call err()
  endif
end subroutine
end module
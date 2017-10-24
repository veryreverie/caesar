! ======================================================================
! I/O operations.
! ======================================================================
module io_module
  use constants_module, only : dp
  use string_module
  implicit none
  
  private
  
  integer      :: TERMINAL_WIDTH = 79
  type(String) :: HOME
  type(String) :: CWD
  type(String) :: OLD_PATH
  type(String) :: PYTHON_PATH
  
  ! Public types.
  public :: CommandLineFlag
  
  ! File operations.
  public :: file_exists ! checks if a file exists
  
  ! Other IO operations.
  public :: set_global_io_variables ! Sets io global variables.
  public :: system_call             ! Makes a system call.
  public :: get_flag                ! Reads a flag from the command line.
  public :: read_line_from_user     ! Reads a line from the terminal.
  public :: print_line              ! write(*,'(a)')
  public :: err                     ! Aborts with a stacktrace.
  public :: format_path             ! Converts any path into an absolute path.
  public :: execute_old_code        ! Runs one of the old caesar codes.
  public :: execute_python          ! Runs one of the python scripts.
  
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
    module procedure print_line_character
    module procedure print_line_String
    module procedure print_line_Stringable
    module procedure print_line_Printable
    
    module procedure print_line_integer
    module procedure print_line_real
    module procedure print_line_logical
    module procedure print_line_complex
    
    module procedure print_line_integers
    module procedure print_line_reals
    module procedure print_line_logicals
    module procedure print_line_complexes
  end interface
  
  interface err
    module procedure err_none
    module procedure err_allocate_flag
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
  
  ! C tput cols interface.
  interface
    function get_terminal_width_c(width) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int), intent(out) :: width
      logical(kind=c_bool)             :: success
    end function
  end interface

contains

! ----------------------------------------------------------------------
! Checks if a file exists
! ----------------------------------------------------------------------
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

! ----------------------------------------------------------------------
! Non-file IO operations.
! ----------------------------------------------------------------------

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
  use iso_fortran_env, only : input_unit
  implicit none
  
  type(String) :: line
  
  character(1000) :: char_line
  
  read(input_unit,'(a)') char_line
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

function get_terminal_width() result(output)
  implicit none
  
  integer :: output
  
  integer :: width
  logical :: success
  
  success = get_terminal_width_c(width)
  if (.not. success) then
    call print_line( 'Failed to get terminal width. &
                     &Reverting to default of 79 characters.')
    output = TERMINAL_WIDTH
  else
    output = width
  endif
end function

subroutine set_global_io_variables()
  implicit none
  
  type(String) :: exe_location
  type(String) :: caesar_dir
  
  TERMINAL_WIDTH = get_terminal_width()
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
  PYTHON_PATH = caesar_dir//'/python'
end subroutine

! ----------------------------------------------------------------------
! Subroutines to print lines, as write(), but with error checking
!    and formatting.
! ----------------------------------------------------------------------
recursive subroutine print_line_character(line)
  use iso_fortran_env, only : output_unit
  implicit none
  
  character(*), intent(in) :: line
  
  integer :: ierr
  integer :: i
  logical :: success
  
  if (len(line) <= TERMINAL_WIDTH) then
    ! The string can be printed all on one line.
    write(output_unit,'(a)',iostat=ierr) line
  else
    success = .false.
    ! Attempt to break the string at a space, so that it fits on the terminal.
    do i=TERMINAL_WIDTH,2,-1
      if (line(i:i)==' ') then
        write(output_unit,'(a)',iostat=ierr) line(:i-1)
        call print_line(line(i+1:))
        success = .true.
        exit
      endif
    enddo
    
    ! Attempt to break the string at the first available space.
    ! Some wrapping will happen, but this can't be avoided.
    if (.not. success) then
      do i=TERMINAL_WIDTH+1,len(line)-1
        if (line(i:i)==' ') then
          write(output_unit,'(a)',iostat=ierr) line(:i-1)
          call print_line(line(i+1:))
          success = .true.
          exit
        endif
      enddo
    endif
    
    ! The line can't be split. Everything is written.
    if (.not. success) then
      write(output_unit,'(a)',iostat=ierr) line
    endif
  endif
  
  if (ierr /= 0) then
    write(output_unit,'(a)') 'Error in print_line.'
    call err()
  endif
  
  flush(output_unit,iostat=ierr)
  
  if (ierr /= 0) then
    write(output_unit,'(a)') 'Error in print_line.'
    call err()
  endif
end subroutine

subroutine print_line_String(line)
  implicit none
  
  type(String), intent(in) :: line
  
  call print_line(char(line))
end subroutine

subroutine print_line_Stringable(this)
  use stringable_module
  implicit none
  
  class(Stringable), intent(in) :: this
  
  call print_line(str(this))
end subroutine

subroutine print_line_Printable(this)
  use printable_module
  implicit none
  
  class(Printable), intent(in) :: this
  
  type(String), allocatable :: lines(:)
  integer                   :: i
  
  lines = this%str()
  do i=1,size(lines)
    call print_line(lines(i))
  enddo
end subroutine

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
  
  call print_line(''//this)
end subroutine

subroutine print_line_reals(this)
  implicit none
  
  real(dp), intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_logicals(this)
  implicit none
  
  logical, intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_complexes(this)
  implicit none
  
  complex(dp), intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

! ----------------------------------------------------------------------
! Aborts with a stacktrace.
! ----------------------------------------------------------------------
! Always aborts.
subroutine err_none()
  use compiler_specific_module
  implicit none
  
  call err_implementation()
end subroutine

! Aborts if integer input /= 0.
! Designed for use with allocate 'stat=ierr' flags.
subroutine err_allocate_flag(this)
  implicit none
  
  integer, intent(in) :: this
  
  if (this/=0) then
    call print_line('Allocation error.')
    call err()
  endif
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
subroutine execute_python(wd, filename)
  implicit none
  
  type(String), intent(in) :: wd
  type(String), intent(in) :: filename
  
  integer :: result_code
  
  result_code = system_call( 'cd '//wd//'; &
                            &python3 '//PYTHON_PATH//'/'//filename)
  
  if (result_code/=0) then
    call print_line('Error: '//filename//' failed.')
    call err()
  endif
end subroutine
end module

! ======================================================================
! I/O operations.
! ======================================================================
module io_module
  use iso_fortran_env,  only : output_unit
  use constants_module, only : dp
  use error_module, only : error_module_ERROR => ERROR, &
                         & error_module_CODE_ERROR => CODE_ERROR, &
                         & error_module_WARNING => WARNING
  use string_module
  implicit none
  
  private
  
  ! The terminal escape character.
  character(1), parameter :: ESC = achar(27)
  
  ! Public strings for ease of printing.
  ! Re-export error strings from error_module.
  character(14), parameter, public :: ERROR = error_module_ERROR
  character(19), parameter, public :: CODE_ERROR = error_module_CODE_ERROR
  character(16), parameter, public :: WARNING = error_module_WARNING
  
  ! I/O settings, specifying various input/output properties.
  ! Set by set_io_settings.
  integer      :: TERMINAL_WIDTH = 79
  type(String) :: HOME
  type(String) :: CWD
  type(String) :: OLD_PATH
  type(String) :: PYTHON_PATH
  integer      :: OUTPUT_FILE_UNIT = output_unit
  
  ! Public types.
  public :: CommandLineFlag
  
  ! Conversions from String.
  public :: lgcl
  public :: int
  public :: dble
  public :: cmplx
  
  ! File operations.
  public :: file_exists ! checks if a file exists
  
  ! Other IO operations.
  public :: set_io_settings     ! Sets I/O settings.
  public :: set_output_unit     ! Sets output file unit.
  public :: unset_output_unit   ! Reverts output file unit to stdout.
  public :: system_call         ! Makes a system call.
  public :: get_flag            ! Reads a flag from the command line.
  public :: read_line_from_user ! Reads a line from the terminal.
  public :: err                 ! Aborts with a stacktrace.
  public :: print_line          ! write(*,'(a)')
  public :: format_path         ! Converts any path into an absolute path.
  public :: execute_old_code    ! Runs one of the old caesar codes.
  public :: execute_python      ! Runs one of the python scripts.
  public :: colour              ! Adds terminal colours to a string.
  public :: left_pad            ! Left pads an integer.
  public :: spaces              ! spaces(n) returns n spaces.
  
  ! Command line flag and argument.
  type :: CommandLineFlag
    character(1) :: flag
    type(String) :: argument
  end type
  
  interface lgcl
    module procedure lgcl_character
    module procedure lgcl_String
  end interface
  
  interface int
    module procedure int_character
    module procedure int_String
  end interface
  
  interface dble
    module procedure dble_character
    module procedure dble_String
  end interface
  
  interface cmplx
    module procedure cmplx_character
    module procedure cmplx_String
  end interface
  
  interface file_exists
    module procedure file_exists_character
    module procedure file_exists_string
  end interface
  
  interface err
    module procedure err_none
    module procedure err_allocate_flag
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
  
  interface format_path
    module procedure format_path_character
    module procedure format_path_String
  end interface
  
  interface colour
    module procedure colour_character_character
    module procedure colour_String_character
    module procedure colour_character_String
    module procedure colour_String_String
  end interface
  
  interface left_pad
    module procedure left_pad_character
    module procedure left_pad_String
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
! Conversion from String
! ----------------------------------------------------------------------
! Convert a character or String to logical.
function lgcl_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  logical                  :: output
  
  integer :: ierr
  
  read(this,*,iostat=ierr) output
  if (ierr/=0) then
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &a logical.')
    call err()
  endif
end function

impure elemental function lgcl_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  logical                  :: output
  
  output = lgcl(char(this))
end function

! Convert a character or String to integer.
function int_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer                  :: output
  
  integer :: ierr
  
  read(this,*,iostat=ierr) output
  if (ierr/=0) then
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &an integer.')
    call err()
  endif
end function

impure elemental function int_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  output = int(char(this))
end function

! Convert a character or String to real(dp).
! Can also convert fractions as strings to real(dp), e.g. '4/5' -> 0.8_dp.
recursive function dble_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  real(dp)                 :: output
  
  type(String), allocatable :: split_string(:)
  real(dp)                  :: reals(2)
  
  integer :: ierr
  
  split_string = split(this,'/')
  
  if (size(split_string)==1) then
    ! The string does not contain a '/': treat it as a real.
    read(this,*,iostat=ierr) output
    if (ierr/=0) then
      call print_line(ERROR//': unable to convert the string "'//this//'" to &
         &a real number.')
      call err()
    endif
  elseif (size(split_string)==2) then
    ! The string contains a '/': treat it as a fraction.
    output = dble(char(split_string(1)))/dble(char(split_string(2)))
  else
    ! The string contains multiple '/'s: this is an error.
    call print_line(ERROR//': unable to convert the string "'//this//'" to &
       &a real number.')
    call err()
  endif
end function

impure elemental function dble_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  real(dp)                 :: output
  
  output = dble(char(this))
end function

! Convert a character or String to complex(dp).
function cmplx_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  complex(dp)              :: output
  
  integer :: i
  logical :: split_allowed
  
  split_allowed = .false.
  do i=len(this)-1,1,-1
    if (this(i:i)=='E') then
      split_allowed = .true.
    elseif (split_allowed .and. this(i:i)=='+') then
      output = cmplx(  dble(this(1:i-1)),           &
                    &  dble(this(i+1:len(this)-1)), &
                    &  dp)
      exit
    elseif (split_allowed .and. this(i:i)=='-') then
      output = cmplx(  dble(this(1:i-1)),           &
                    & -dble(this(i+1:len(this)-1)), &
                    &  dp)
      exit
    endif
  enddo
end function

impure elemental function cmplx_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  complex(dp)              :: output
  
  output = cmplx(char(this))
end function

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

subroutine set_io_settings()
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
! Set or unset OUTPUT_FILE_UNIT as a file unit.
! ----------------------------------------------------------------------
subroutine set_output_unit(file_unit)
  implicit none
  
  integer, intent(in) :: file_unit
  
  if (OUTPUT_FILE_UNIT/=output_unit) then
    call print_line('Code Error: attempted to redirect stdout when it was &
       &already redirected.')
    call err()
  endif
  
  OUTPUT_FILE_UNIT = file_unit
end subroutine

subroutine unset_output_unit()
  implicit none
  
  OUTPUT_FILE_UNIT = output_unit
end subroutine

! ----------------------------------------------------------------------
! Aborts with a stacktrace.
! ----------------------------------------------------------------------
! Always aborts.
subroutine err_none()
  use error_module, only : abort_with_stacktrace
  implicit none
  
  call abort_with_stacktrace()
end subroutine

! Aborts if integer input /= 0.
! Designed for use with allocate 'stat=ierr' flags.
subroutine err_allocate_flag(this)
  implicit none
  
  integer, intent(in) :: this
  
  if (this/=0) then
    write(*,*) 'Allocation error.'
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! Subroutines to print lines, as write(), but with error checking
!    and formatting.
! ----------------------------------------------------------------------
subroutine print_line_character(line,indent)
  implicit none
  
  character(*), intent(in)           :: line
  integer,      intent(in), optional :: indent
  
  integer              :: line_length
  integer              :: no_lines
  integer, allocatable :: line_starts(:)
  integer, allocatable :: line_ends(:)
  
  integer                   :: last_space
  logical                   :: reading_special
  integer                   :: current_terminal_width
  character(:), allocatable :: indent_spaces
  character(:), allocatable :: overhang_spaces
  
  integer :: i,ierr,ialloc
  
  if (present(indent)) then
    indent_spaces = spaces(indent)
  else
    indent_spaces = spaces(0)
  endif
  overhang_spaces = spaces(3)
  
  allocate( line_starts(max(len(line),1)), &
          & line_ends(max(len(line),1)),   &
          & stat=ialloc); call err(ialloc)
  
  current_terminal_width = TERMINAL_WIDTH - len(indent_spaces)
  no_lines = 1
  line_starts(1) = 1
  last_space = 0
  
  ! If not writing to the terminal, ignore line-breaking.
  if (OUTPUT_FILE_UNIT/=output_unit) then
    line_ends(1) = len(line)
  
  ! If writing to the terminal, identify line breaks.
  else
    reading_special = .false.
    line_length = 0
    do i=1,len(line)
      ! Special characters, which don't display on the terminal.
      if (line(i:i)==ESC) then
        reading_special = .true.
      elseif (reading_special) then
        if (line(i:i)=='m') then
          reading_special = .false.
        endif
      
      ! Normal characters, which do display on the terminal.
      else
        line_length = line_length+1
        if (line_length>current_terminal_width) then
          if (last_space/=0) then
            line_ends(no_lines) = last_space
            line_starts(no_lines+1) = last_space+1
          else
            line_ends(no_lines) = i-1
            line_starts(no_lines+1) = i
          endif
          line_length = line_length &
                    & - (line_ends(no_lines)-line_starts(no_lines)+1)
          no_lines = no_lines+1
          last_space = 0
          current_terminal_width = TERMINAL_WIDTH     &
                               & - len(indent_spaces) &
                               & - len(overhang_spaces)
        elseif (line(i:i)==' ') then
          last_space = i
        endif
      endif
    enddo
  endif
  
  line_ends(no_lines) = len(line)
  
  ! Write lines.
  do i=1,no_lines
    if (i==1) then
      write(OUTPUT_FILE_UNIT,'(a)',iostat=ierr) &
         & indent_spaces//line(line_starts(i):line_ends(i))
    else
      write(OUTPUT_FILE_UNIT,'(a)',iostat=ierr) &
         & indent_spaces//overhang_spaces//line(line_starts(i):line_ends(i))
    endif
    
    ! Check for errors and flush the write buffer.
    if (ierr /= 0) then
      write(OUTPUT_FILE_UNIT,'(a)') 'Error in print_line.'
      call err()
    endif
    
    flush(OUTPUT_FILE_UNIT,iostat=ierr)
    
    if (ierr /= 0) then
      write(OUTPUT_FILE_UNIT,'(a)') 'Error in print_line.'
      call err()
    endif
  enddo
end subroutine

subroutine print_line_String(line,indent)
  implicit none
  
  type(String), intent(in)           :: line
  integer,      intent(in), optional :: indent
  
  if (present(indent)) then
    call print_line(char(line),indent)
  else
    call print_line(char(line))
  endif
end subroutine

subroutine print_line_Stringable(this,indent)
  implicit none
  
  class(Stringable), intent(in)           :: this
  integer,           intent(in), optional :: indent
  
  if (present(indent)) then
    call print_line(str(this),indent)
  else
    call print_line(str(this))
  endif
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

function spaces(no_spaces) result(output)
  implicit none
  
  integer, intent(in)       :: no_spaces
  character(:), allocatable :: output
  
  integer :: i,ialloc
  
  allocate(character(no_spaces) :: output, stat=ialloc); call err(ialloc)
  
  do i=1,no_spaces
    output(i:i) = ' '
  enddo
end function

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

! ----------------------------------------------------------------------
! Adds terminal colour to a string.
! Does nothing if the output is going to a file.
! ----------------------------------------------------------------------
! N.B. can't detect "caesar > file", only supresses colour if "caesar -o file".
function colour_character_character(input,colour_name) result(output)
  implicit none
  
  character(*), intent(in) :: input
  character(*), intent(in) :: colour_name
  type(String)             :: output
  
  type(String) :: lower_case_name
  character(2) :: colour_code
  
  ! Supress colour if outputting to file.
  if (OUTPUT_FILE_UNIT/=output_unit) then
    output = input
    return
  endif
  
  lower_case_name = lower_case(colour_name)
  
  if (lower_case_name=='black') then
    colour_code = '30'
  elseif (lower_case_name=='red') then
    colour_code = '31'
  elseif (lower_case_name=='green') then
    colour_code = '32'
  elseif (lower_case_name=='yellow') then
    colour_code = '33'
  elseif (lower_case_name=='blue') then
    colour_code = '34'
  elseif (lower_case_name=='magenta') then
    colour_code = '35'
  elseif (lower_case_name=='cyan') then
    colour_code = '36'
  elseif (lower_case_name=='light gray') then
    colour_code = '37'
  elseif (lower_case_name=='dark gray') then
    colour_code = '90'
  elseif (lower_case_name=='light red') then
    colour_code = '91'
  elseif (lower_case_name=='light green') then
    colour_code = '92'
  elseif (lower_case_name=='light yellow') then
    colour_code = '93'
  elseif (lower_case_name=='light blue') then
    colour_code = '94'
  elseif (lower_case_name=='light magenta') then
    colour_code = '95'
  elseif (lower_case_name=='light cyan') then
    colour_code = '96'
  elseif (lower_case_name=='white') then
    colour_code = '97'
  else
    call print_line('Error: '//colour_name//' is not an accepted colour.')
    call err()
  endif
  
  output = ESC//'['//colour_code//'m'//input//ESC//'[0m'
end function

function colour_String_character(input,colour_name) result(output)
  implicit none
  
  character(*), intent(in) :: input
  type(String), intent(in) :: colour_name
  type(String)             :: output
  
  output = colour(input,char(colour_name))
end function

function colour_character_String(input,colour_name) result(output)
  implicit none
  
  type(String), intent(in) :: input
  character(*), intent(in) :: colour_name
  type(String)             :: output
  
  output = colour(char(input),colour_name)
end function

function colour_String_String(input,colour_name) result(output)
  implicit none
  
  type(String), intent(in) :: input
  type(String), intent(in) :: colour_name
  type(String)             :: output
  
  output = colour(char(input),char(colour_name))
end function

! ----------------------------------------------------------------------
! Left pads an integer, by default with zeroes.
! ----------------------------------------------------------------------

! Left pads to match the length of a given string.
! e.g. left_pad(5,'example') returns '0000005'.
! Throws an error if the input is negative, or longer than the matching string.
function left_pad_character(input,match,pad_character) result(output)
  implicit none
  
  integer,      intent(in)           :: input
  character(*), intent(in)           :: match
  character(1), intent(in), optional :: pad_character
  type(String)                       :: output
  
  character(1) :: pad
  integer      :: i
  
  if (present(pad_character)) then
    pad = pad_character
  else
    pad = '0'
  endif
  
  if (input<0) then
    call err()
  endif
  
  output = str(input)
  
  if (len(output)>len(match)) then
    call err()
  endif
  
  do i=1,len(match)-len(output)
    output = pad//output
  enddo
end function

! As above.
function left_pad_String(input,match,pad_character) result(output)
  implicit none
  
  integer,      intent(in)           :: input
  type(String), intent(in)           :: match
  character(1), intent(in), optional :: pad_character
  type(String)                       :: output
  
  if (present(pad_character)) then
    output = left_pad(input, char(match), pad_character)
  else
    output = left_pad(input, char(match))
  endif
end function
end module

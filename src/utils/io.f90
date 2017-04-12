module io_module
  use constants_module, only : dp
  use string_module
  implicit none
  
  private
  
  integer :: TERMINAL_WIDTH = 79
  
  ! Command line flag and argument.
  public CommandLineFlag
  
  ! File operations.
  public :: open_read_file   ! open a file for reading
  public :: open_write_file  ! open a file for writing
  public :: open_append_file ! open a file for appending to
  public :: file_exists      ! checks if a file exists
  public :: count_lines      ! counts the number of lines in a file
  public :: read_lines       ! reads a file into a String(:) array
  
  ! Other IO operations.
  public :: system_call           ! Makes a system call.
  public :: get_current_directory ! Gets the current directory.
  public :: get_flag              ! Reads a flag from the command line.
  public :: read_line_from_user   ! Reads a line from the terminal.
  public :: update_terminal_width ! Gets the terminal width.
  public :: print_line            ! write(*,'(a)')
  public :: err                   ! Aborts with a stacktrace.
  
  type CommandLineFlag
    character(1) :: flag
    type(String) :: argument
  end type
  
  interface open_read_file
    module procedure open_read_file_character
    module procedure open_read_file_string
  end interface
  
  interface open_write_file
    module procedure open_write_file_character
    module procedure open_write_file_string
  end interface
  
  interface open_append_file
    module procedure open_append_file_character
    module procedure open_append_file_string
  end interface
  
  interface file_exists
    module procedure file_exists_character
    module procedure file_exists_string
  end interface
  
  interface count_lines
    module procedure count_lines_character
    module procedure count_lines_string
  end interface
  
  interface read_lines
    module procedure read_lines_character
    module procedure read_lines_String
  end interface
  
  interface print_line
    module procedure print_line_character
    module procedure print_line_file_character
    module procedure print_line_String
    module procedure print_line_file_String
    
    module procedure print_line_integer
    module procedure print_line_file_integer
    module procedure print_line_real
    module procedure print_line_file_real
    module procedure print_line_logical
    module procedure print_line_file_logical
    
    module procedure print_line_integers
    module procedure print_line_file_integers
    module procedure print_line_reals
    module procedure print_line_file_reals
    module procedure print_line_logicals
    module procedure print_line_file_logicals
  end interface
  
  interface err
    module procedure err_none
    module procedure err_logical
    module procedure err_integer
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
  
  ! C pwd call interface.
  interface
    function pwd_c(cwd_size, cwd) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      implicit none
      
      integer(kind=c_int),    intent(in)  :: cwd_size
      character(kind=c_char), intent(out) :: cwd(*)
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

! open a file with a specified mode, and return the unit it is opened in
function open_file(filename,status,action,access) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  character(*), intent(in) :: status
  character(*), intent(in) :: action
  character(*), intent(in) :: access
  integer                  :: unit_num
  
  integer :: iostat
  logical :: opened
  
  unit_num = 0
  opened=.true.
  
  do while (opened) ! loop until unused unit is found
    unit_num = unit_num + 1
    if (unit_num==5 .or. unit_num==6) cycle ! ignore 0,5,6 std units
    if (unit_num>=100 .and. unit_num<=102) cycle ! ignore 100,101,102 std units
    
    inquire(unit=unit_num,opened=opened,iostat=iostat) ! check unit
    if (iostat/=0) opened=.true. ! ignore unit if error
  enddo
  
  open(unit=unit_num,file=filename,status=status,action=action,&
    &access=access,iostat=iostat)
  if (iostat /= 0) then
    call print_line('Error opening '//filename//' file.')
    call err()
  endif
end function

! open a file for reading, and return the unit it is opened in
function open_read_file_character(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'old', 'read', 'sequential')
end function

function open_read_file_string(filename) result(unit_num)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_read_file(char(filename))
end function

! open a file for writing, and return the unit it is opened in
function open_write_file_character(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'sequential')
end function

function open_write_file_string(filename) result(unit_num)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_write_file(char(filename))
end function

! open a file for appending, and return the unit it is opened in
function open_append_file_character(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'append')
end function

function open_append_file_string(filename) result(unit_num)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_append_file(char(filename))
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
! Gets the number of lines remaining in a file
! ----------------------------------------------------------------------
function count_lines_character(filename) result(output)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: output
  
  integer      :: file_unit
  integer      :: iostat
  character(1) :: line
  
  file_unit = open_read_file(filename)
  output = 0
  iostat = 0
  do while (iostat==0)
    read(file_unit, '(a)', iostat=iostat) line
    if (iostat==0) then
      output = output+1
    elseif (iostat>0) then
      call print_line('Error counting lines of '//filename)
      call err()
    endif
  enddo
  close(file_unit)
end function

function count_lines_String(filename) result(output)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: output
  
  output = count_lines(char(filename))
end function

! ----------------------------------------------------------------------
! Reads the file into a String(:) array
! ----------------------------------------------------------------------
function read_lines_character(filename) result(output)
  implicit none
  
  character(*), intent(in)  :: filename
  type(String), allocatable :: output(:)
  
  integer         :: file_length
  integer         :: file_unit
  character(1000) :: line
  integer         :: i
  
  integer :: ierr
  
  file_length = count_lines(filename)
  
  allocate(output(file_length),stat=ierr)
  if (ierr/=0) then
    call print_line('Error allocating output for read_lines')
    call err()
  endif
  
  file_unit = open_read_file(filename)
  do i=1,file_length
    read(file_unit,'(a)',iostat=ierr) line
    if (ierr/=0) then
      call print_line('Error reading from '//filename)
      call err()
    endif
    output(i) = trim(line)
  enddo
  close(file_unit)
end function

function read_lines_String(filename) result(output)
  implicit none
  
  type(String), intent(in)  :: filename
  type(String), allocatable :: output(:)
  
  output = read_lines(char(filename))
end function

! ----------------------------------------------------------------------
! Non-file IO operations.
! ----------------------------------------------------------------------

! Makes a system call via system.c.
subroutine system_call(this)
  implicit none
  
  type(String), intent(in) :: this
  
  integer :: result_code
  
  result_code = system_c(char(this)//char(0))
end subroutine

! Gets the current working directory via system.c.
function get_current_directory() result(output)
  implicit none
  
  type(String) :: output
  
  integer, parameter  :: cwd_size = 1024
  
  character(cwd_size) :: current_dir
  logical             :: success
  integer             :: i
  
  success = pwd_c(cwd_size, current_dir)
  
  if (.not. success) then
    call print_line('pwd failed.')
    call err()
  endif
  
  output = ''
  do i=1,cwd_size
    if (current_dir(i:i)==char(0)) then
      output = current_dir(:i-1)
      exit
    endif
  enddo
  
  if (output=='') then
    call print_line('pwd string not nul-terminated.'//current_dir)
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Reads flags from the command line via system.c.
! ----------------------------------------------------------------------
! Returns none-flag arguments as flag=' '
! Aborts if an unexpected flag is passed.
! Aborts if no argument is passed where one is required.
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
  argvs = char('caesar'//char(0)//join(args,char(0))//char(0))
  
  ! Convert flags into C-friendly format.
  options = '-:'//flags_without_arguments
  flags_char = char(flags_with_arguments)
  do i=1,len(flags_with_arguments)
    options = options//flags_char(i:i)//':'
  enddo
  options = options//char(0)
  
  ! Call getopt.
  success = get_flag_c( size(args)+1,  &
                      & argvs,         &
                      & char(options), &
                      & flag,          &
                      & argument)
  
  ! Return result.
  if (.not. success) then
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
    
    if (output%flag=='?') then
      call print_line('Error: unexpected flag "'//output%argument//'".')
      stop
    elseif (output%flag==':') then
      call print_line('Error: no option given for flag "'// &
         & output%argument//'".')
      stop
    endif
  endif
end function

function read_line_from_user() result(line)
  implicit none
  
  type(String) :: line
  
  character(1000) :: char_line
  
  read(*,'(a)') char_line
  line = trim(char_line)
end function

! ----------------------------------------------------------------------
! Checks the width of the terminal, and sets TERMINAL_WIDTH.
! ----------------------------------------------------------------------
subroutine update_terminal_width()
  implicit none
  
  logical :: success
  integer :: width
  
  success = get_terminal_width_c(width)
  if (.not. success) then
    call print_line( 'Failed to get terminal width. &
                     &Reverting to default of 79 characters.')
  else
    TERMINAL_WIDTH = width
  endif
end subroutine

! ----------------------------------------------------------------------
! Subroutines to print lines, as write(), but with error checking
!    and formatting.
! ----------------------------------------------------------------------
recursive subroutine print_line_character(line)
  implicit none
  
  character(*), intent(in) :: line
  
  integer :: ierr
  integer :: i
  logical :: success
  
  if (len(line) <= TERMINAL_WIDTH) then
    ! The string can be printed all on one line.
    write(*,'(a)',iostat=ierr) line
  else
    success = .false.
    ! Attempt to break the string at a space, so that it fits on the terminal.
    do i=TERMINAL_WIDTH,2,-1
      if (line(i:i)==' ') then
        write(*,'(a)',iostat=ierr) line(:i-1)
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
          write(*,'(a)',iostat=ierr) line(:i-1)
          call print_line(line(i+1:))
          success = .true.
          exit
        endif
      enddo
    endif
    
    ! The line can't be split. Everything is written.
    if (.not. success) then
      write(*,'(a)',iostat=ierr) line
    endif
  endif
  
  if (ierr /= 0) then
    write(*,*) 'Error in print_line.'
    call err()
  endif
end subroutine

subroutine print_line_file_character(file_unit,line)
  implicit none
  
  integer,      intent(in) :: file_unit
  character(*), intent(in) :: line
  
  integer :: ierr
  
  write(file_unit,'(a)',iostat=ierr) line
  
  if (ierr /= 0) then
    write(*,*) 'Error in print_line.'
    call err()
  endif
end subroutine

subroutine print_line_String(line)
  implicit none
  
  type(String), intent(in) :: line
  
  call print_line(char(line))
end subroutine

subroutine print_line_file_String(file_unit,line)
  implicit none
  
  integer,      intent(in) :: file_unit
  type(String), intent(in) :: line
  
  call print_line(file_unit,char(line))
end subroutine

subroutine print_line_integer(this)
  implicit none
  
  integer, intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_integer(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  integer, intent(in) :: this
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_real(this)
  implicit none
  
  real(dp), intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_real(file_unit,this)
  implicit none
  
  integer,  intent(in) :: file_unit
  real(dp), intent(in) :: this
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_logical(this)
  implicit none
  
  logical, intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_logical(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  logical, intent(in) :: this
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_integers(this)
  implicit none
  
  integer, intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_integers(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  integer, intent(in) :: this(:)
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_reals(this)
  implicit none
  
  real(dp), intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_reals(file_unit,this)
  implicit none
  
  integer,  intent(in) :: file_unit
  real(dp), intent(in) :: this(:)
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_logicals(this)
  implicit none
  
  logical, intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_logicals(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  logical, intent(in) :: this(:)
  
  call print_line(file_unit,''//this)
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

! Aborts if logical input is .false.
subroutine err_logical(this)
  implicit none
  
  logical, intent(in) :: this
  
  if (.not. this) then
    call err()
  endif
end subroutine

! Aborts if integer input /= 0.
! Designed for use with allocate 'stat=ialloc' flags.
subroutine err_integer(this)
  implicit none
  
  integer, intent(in) :: this
  
  if (this/=0) then
    call print_line('Allocation error.')
    call err()
  endif
end subroutine
end module

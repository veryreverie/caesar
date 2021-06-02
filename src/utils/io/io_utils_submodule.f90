submodule (caesar_io_utils_module) caesar_io_utils_submodule
  use caesar_home_directory_module
  use caesar_current_working_directory_module
  use caesar_executable_locations_module
contains

module procedure file_exists_character
  inquire(file=char(format_path(filename)), exist=output)
end procedure

module procedure file_exists_string
  output = file_exists(char(filename))
end procedure

module procedure command_line_args
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
end procedure

module procedure mkdir
  type(String) :: formatted_dirname
  
  integer :: result_code
  
  formatted_dirname = format_path(dirname)
  
  result_code = mkdir_c(char(formatted_dirname)//char(0))
  
  if (result_code/=0) then
    call print_line(ERROR//': failed to make directory: '//formatted_dirname)
    call print_line('mkdir return code: '//result_code)
    call err()
  endif
end procedure

module procedure system_call
  result_code = system_c(char(this)//char(0))
end procedure

module procedure get_flag
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
      call print_line(ERROR//': unexpected flag "'//output%argument//'".')
      call quit()
    
    ! The flag has no argument.
    elseif (output%flag==':') then
      output%flag = char(output%argument)
      output%argument = ''
    endif
  endif
end procedure

module procedure read_line_from_user
  character(1000) :: char_line
  
  read(INPUT_UNIT,'(a)') char_line
  line = trim(char_line)
end procedure

module procedure format_path_character
  integer :: length
  
  type(String) :: wd
  
  length = len(path)
  
  if (length==0) then
    call print_line(ERROR//': no path provided.')
    call err()
  endif
  
  if (present(working_directory)) then
    wd = working_directory
  else
    wd = current_working_directory()
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
      output = wd
    elseif (path(:1)=='~') then
      ! ~ -> /home/directory
      output = home_directory()
    else
      ! foo -> /current/working/directory/foo
      output = wd//'/'//path(:1)
    endif
  else
    if (path(:1)=='/') then
      ! /foo -> /foo
      output = path(:length)
    elseif (path(:2)=='./') then
      ! ./foo -> /current/working/directory/foo
      output = wd//path(2:length)
    elseif (path(:2)=='~/') then
      ! ~/foo -> /home/directory/foo
      output = home_directory()//path(2:length)
    else
      ! foo -> /current/working/directory/foo
      output = wd//'/'//path(:length)
    endif
  endif
end procedure

module procedure format_path_String
  output = format_path(char(path),working_directory)
end procedure

module procedure execute_python
  type(String) :: command
  integer      :: result_code
  
  command = 'cd '//current_working_directory()//'; &
            &'//python_path//' '//python_scripts_location()//'/'//filename
  
  if (present(python_arguments)) then
    command = command//' '//join(python_arguments)
  endif
  
  result_code = system_call(command)
  
  if (result_code/=0) then
    call print_line(ERROR//': '//filename//' failed.')
    call err()
  endif
end procedure

module procedure call_caesar_character
  type(String) :: command
  integer      :: result_code
  
  command = 'cd '//current_working_directory()//'; ' &
       & // executable_location()//'/caesar'
  if (present(arguments)) then
    command = command//' '//escape_bash_characters(arguments)
  endif
  
  result_code = system_call(command)
  
  if (result_code/=0) then
    call print_line(ERROR//': Calling Caesar from within Caesar failed.')
    call print_line('Caesar was called as:')
    call print_line(command)
    call err()
  endif
end procedure

module procedure call_caesar_String
  call call_caesar(char(arguments))
end procedure

module procedure escape_bash_characters
  integer :: i,j
  
  output = ''
  
  j = 0
  do i=1,len(input)
    if (input(i:i)=='|' .or. input(i:i)=='\') then
      output = output//input(j+1:i-1)//'\'//input(i:i)
      j = i
    endif
  enddo
  output = output//input(j+1:)
end procedure
end submodule

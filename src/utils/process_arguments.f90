! ======================================================================
! Processes arguments from command line, file and interactive input.
! Stops with an error message if unexpected input is given.
! ======================================================================
module process_arguments_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  abstract interface
    subroutine MainSubroutine(arguments)
      use dictionary_module
      implicit none
      
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
contains

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
function process_arguments(mode,args,keywords) result(arguments)
  use dictionary_module
  use help_module
  implicit none
  
  type(String),      intent(in) :: mode
  type(String),      intent(in) :: args(:)
  type(KeywordData), intent(in) :: keywords(:)
  type(Dictionary)              :: arguments
  
  ! Flags.
  type(String)              :: flags_without_arguments
  type(String), allocatable :: long_flags_without_arguments(:)
  
  type(String)              :: flags_with_arguments
  type(String), allocatable :: long_flags_with_arguments(:)
  type(String), allocatable :: default_arguments(:)
  
  type(CommandLineFlag)     :: flag
  
  ! Non-flag inputs.
  integer                   :: no_args
  type(String), allocatable :: arg_keys(:)
  type(String), allocatable :: arg_values(:)
  type(Dictionary)          :: command_line_arguments
  type(Dictionary)          :: file_arguments
  type(String)              :: wd
  type(String)              :: input_filename
  logical                   :: interactive
  
  ! Temporary variables.
  integer                   :: i,j,ialloc
  type(String)              :: temp_string
  logical                   :: appending
  
  ! --------------------------------------------------
  ! Process arguments and flags.
  ! --------------------------------------------------
  
  ! Initialise flags.
  flags_without_arguments = 'i'
  long_flags_without_arguments = [ str('interactive') ]
  if (len(flags_without_arguments)/=size(long_flags_without_arguments)) then
    call err()
  endif
  
  flags_with_arguments = 'dfh'
  long_flags_with_arguments = [ str('working_directory'), &
                              & str('input_file'),        &
                              & str('help')               ]
  if (len(flags_with_arguments)/=size(long_flags_with_arguments)) then
    call err()
  endif
  default_arguments = [ str('.'), str(NOT_SET), str(NOT_SET) ]
  if (len(flags_with_arguments)/=size(default_arguments)) then
    call err()
  endif
  
  ! Set dictionary of all keywords.
  no_args = size(keywords)                     &
        & + size(long_flags_without_arguments) &
        & + size(long_flags_with_arguments)
  allocate( arg_keys(no_args),   &
          & arg_values(no_args), &
          & stat=ialloc); call err(ialloc)
  
  j = 0
  do i=1,size(long_flags_without_arguments)
    j = j+1
    arg_keys(j) = long_flags_without_arguments(i)
    arg_values(j) = NOT_SET
  enddo
  
  do i=1,size(long_flags_with_arguments)
    j = j+1
    arg_keys(j) = long_flags_with_arguments(i)
    arg_values(j) = default_arguments(i)
  enddo
  
  do i=1,size(keywords)
    j = j+1
    arg_keys(j) = keywords(i)%keyword
    arg_values(j) = keywords(i)%default_value
  enddo
  if (j/=no_args) then
    call err()
  endif
  arguments = make_dictionary(arg_keys, arg_values)
  
  ! --------------------------------------------------
  ! Parse command line arguments.
  ! --------------------------------------------------
  no_args = 0
  
  do_args : do
    flag = get_flag(args, flags_without_arguments, flags_with_arguments)
    
    ! No flags remaining.
    if (flag%flag==' ' .and. flag%argument=='') then
      exit
    
    ! The mode.
    elseif (flag%flag==' ' .and. no_args==0) then
      if (flag%argument/=mode) then
        call print_line('Error: Caesar only takes one none-keyword argument. &
           &all other arguments should be preceded by flags or "--" keywords.')
        stop
      endif
    
    ! An argument which is not preceded by '-' or '--'.
    elseif (flag%flag==' ') then
      if (appending) then
        arg_values(no_args) = arg_values(no_args)//' '//flag%argument
      else
        arg_values(no_args) = flag%argument
        appending = .true.
      endif
    
    ! An argument which is preceded by '--', i.e. a long form option.
    elseif (flag%flag=='-') then
      
      ! Check the key is valid.
      if (index(arguments, flag%argument)==0) then
        call print_line('Error: unexpected command-line argument: '// &
           & flag%argument)
        stop
      endif
      
      ! Add the argument to the argument list.
      no_args = no_args+1
      arg_keys(no_args) = flag%argument
      arg_values(no_args) = NO_ARGUMENT
      appending = .false.
    
    ! An argument which is preceded by '-', i.e. a flag.
    else
      ! Flags which do not accept arguments.
      do i=1,len(flags_without_arguments)
        if (flag%flag==slice(flags_without_arguments,i,i)) then
          no_args = no_args+1
          arg_keys(no_args) = long_flags_without_arguments(i)
          arg_values(no_args) = NO_ARGUMENT
          appending = .false.
          cycle do_args
        endif
      enddo
      
      ! Flags which accept arguments.
      do i=1,len(flags_with_arguments)
        if (flag%flag==slice(flags_with_arguments,i,i)) then
          no_args = no_args+1
          arg_keys(no_args) = long_flags_with_arguments(i)
          if (flag%argument=='') then
            arg_values(no_args) = NO_ARGUMENT
            appending = .false.
          else
            arg_values(no_args) = flag%argument
            appending = .true.
          endif
          cycle do_args
        endif
      enddo
      
    endif
  enddo do_args
  
  ! Check for duplicate arguments.
  do i=1,no_args
    do j=i+1,no_args
      if (arg_keys(i)==arg_keys(j)) then
        call print_line('Error: argument '//arg_keys(i)//' specified more &
           &than once on the command line. This may be a result of also &
           &specifying the equivalent flag.')
        stop
      endif
    enddo
  enddo
  
  ! Convert arguments into a dictionary.
  command_line_arguments = make_dictionary( arg_keys(:no_args), &
                                          & arg_values(:no_args))
  
  ! Check if interactive mode requested.
  interactive = item(command_line_arguments, 'interactive') /= NOT_SET
  
  ! --------------------------------------------------
  ! Process input file, if requested.
  ! --------------------------------------------------
  input_filename = item(command_line_arguments, 'input_file')
  if (input_filename==NO_ARGUMENT) then
    call print_line('Error: no argument given with --input_file.')
    stop
  endif
  
  if (interactive) then
    call print_line('')
    call print_line('Read additional settings from a file? (y/n)')
    temp_string = read_line_from_user()
    
    if (temp_string=='y') then
      call print_line('')
      if (input_filename/=NOT_SET) then
        call print_line('Input file is currently set to: '//input_filename)
        call print_line('Press <Enter> to accept current setting.')
      endif
      call print_line('Where is the input file?')
      temp_string = read_line_from_user()
      if (temp_string/='') then
        input_filename = temp_string
      endif
    else
      input_filename = NOT_SET
    endif
  endif
  
  if (input_filename/=NOT_SET) then
    file_arguments = read_dictionary_file(input_filename)
    
    ! Check all file arguments are expected.
    do i=1,size(file_arguments)
      if (index(arguments, file_arguments%keys(i))==0) then
        call print_line('Error: unexpected keyword '//file_arguments%keys(i) &
           & //' in file '//input_filename//'.')
        stop
      elseif ( file_arguments%keys(i)=='input_filename' .or. &
             & file_arguments%keys(i)=='help'           .or. &
             & file_arguments%keys(i)=='interactive') then
        file_arguments%values(i) = NOT_SET
      endif
    enddo
  endif
  
  ! --------------------------------------------------
  ! Combine file and command line input arguments.
  ! --------------------------------------------------
  ! N.B. Command line arguments overwrite file arguments,
  !    and both overwrite default arguments.
  do i=1,size(arguments)
    temp_string = item(command_line_arguments, arguments%keys(i))
    if (temp_string/=NOT_SET) then
      arguments%values(i) = temp_string
    else if (input_filename/=NOT_SET) then
      temp_string = item(file_arguments, arguments%keys(i))
      if (temp_string/=NOT_SET) then
        arguments%values(i) = temp_string
      endif
    endif
    
    if (.not. interactive) then
      if (arguments%values(i)==NOT_SET) then
        cycle
      elseif ( arguments%keys(i)=='interactive' .and. &
             & arguments%values(i)/=NO_ARGUMENT) then
        call print_line('Error: an argument was given for keyword '// &
           & arguments%keys(i)//', which does not accept arguments.')
        stop
      elseif ( arguments%keys(i)/='interactive' .and. &
             & arguments%values(i)==NO_ARGUMENT) then
        call print_line('Error: no argument given for keyword '// &
           & arguments%keys(i)//', which requires an argument.')
        stop
      endif
    endif
  enddo
  
  ! Set working directory.
  j = index(arguments, 'working_directory')
  wd = arguments%values(j)
  
  if (interactive) then
    call print_line('')
    call print_line('Where is the working directory?')
    call print_line('Working directory is currently set to: '//wd)
    call print_line('Press <Enter> to accept current setting.')
    temp_string = read_line_from_user()
    if (temp_string/='') then
      wd = temp_string
    endif
  endif
  
  wd = format_path(wd)
  arguments%values(j) = wd
  
  ! --------------------------------------------------
  ! Get interactive arguments, if requested.
  ! --------------------------------------------------
  if (interactive) then
    do i=1,size(keywords)
      j = size(arguments)-size(keywords)+i
      
      if (keywords(i)%keyword/=arguments%keys(j)) then
        call err()
      endif
      
      call print_line('')
      do
        call print_line('Please give a value for the keyword: '// &
           & keywords(i)%keyword//'.')
        call print_line(keywords(i)%helptext)
        if (arguments%values(j)/=NO_ARGUMENT) then
          call print_line(keywords(i)%keyword//' is currently set to: '//&
             & arguments%values(j))
          call print_line('Press <Enter> to accept current setting.')
        endif
        temp_string = read_line_from_user()
        if (temp_string/='') then
          arguments%values(j) = temp_string
        endif
        
        if (arguments%values(j)/=NO_ARGUMENT) then
          exit
        endif
      enddo
      
    enddo
  endif
  
  ! --------------------------------------------------
  ! Convert all file/folder paths to absolute format.
  ! --------------------------------------------------
  do i=1,size(keywords)
    j = size(arguments)-size(keywords)+i
    if (keywords(i)%is_path) then
      if (arguments%values(j) == NOT_SET) then
        cycle
      elseif (arguments%values(j) == NO_ARGUMENT) then
        ! This should have been dealt with above.
        call err()
      else
        arguments%values(j) = format_path(arguments%values(j))
      endif
    endif
  enddo
end function

end module

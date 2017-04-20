program caesar
  ! use utility modules
  use string_module
  use io_module
  use dictionary_module
  use utils_module, only : command_line_args, format_path
  use help_module
  
  ! use harmonic modules
  use setup_harmonic_module
  use run_harmonic_module
  use lte_harmonic_module
  
  ! use quadratic modules
  use setup_quadratic_module
  use run_quadratic_module
  use anharmonic_module
  use bs_quadratic_module
  
  ! use testing modules
  use test_copy_harmonic_module
  use test_lte_module
  use test_copy_quadratic_module
  
  ! use misc modules
  use calculate_gap_module
  use hartree_to_eV_module
  
  implicit none
  
  ! Working directory.
  type(String) :: cwd
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  
  type(String)              :: flags_without_arguments
  type(String), allocatable :: long_flags_without_arguments(:)
  
  type(String)              :: flags_with_arguments
  type(String), allocatable :: long_flags_with_arguments(:)
  
  type(CommandLineFlag)     :: flag
  type(String)              :: mode
  
  ! Input arguments.
  integer                   :: no_args
  type(String), allocatable :: arg_keys(:)
  type(String), allocatable :: arg_values(:)
  type(Dictionary)          :: arguments
  type(String), allocatable :: input_file(:)
  
  ! Keywords which the program accepts.
  type(KeywordData), allocatable :: keywords(:)
  
  ! Inputs.
  type(String) :: wd
  type(String) :: input_filename
  logical      :: interactive
  
  ! Temporary variables.
  integer                   :: i,j,ialloc
  type(String)              :: temp_string
  type(String)              :: temp_string_2
  type(String), allocatable :: line(:)
  
  ! Set terminal width for formatting purposes.
  call update_terminal_width()
  
  ! Get current directory.
  cwd = get_current_directory()
  
  ! Read in command line arguments.
  args = command_line_args()
  if (size(args) < 1) then
    call print_line('No arguments given. For help, call caesar -h')
    stop
  endif
  
  ! Allocate space for input arguments.
  no_args = 0
  allocate( arg_keys(size(args)), &
          & arg_values(size(args)), &
          & stat=ialloc); call err(ialloc)
  
  ! Set default mode.
  mode = ''
  
  ! Set up space for flag handling.
  flags_without_arguments = 'hi'
  long_flags_without_arguments = [ str('help'), str('interactive') ]
  
  flags_with_arguments = 'df'
  long_flags_with_arguments = [ str('working_directory'), str('input_file') ]
  
  ! Parse command line arguments.
  do_args : do
    flag = get_flag(args, flags_without_arguments, flags_with_arguments)
    
    ! No flags remaining.
    if (flag%flag==' ' .and. flag%argument=='') then
      exit
    
    ! The mode.
    elseif (flag%flag==' ' .and. no_args==0) then
      if (mode/='') then
        call print_line('Error: more than one mode specified.')
        stop
      endif
      
      mode = flag%argument
    
    ! An argument which is not preceded by '-' or '--'.
    elseif (flag%flag==' ') then
      if (arg_values(no_args)=='<no_argument>') then
        arg_values(no_args) = flag%argument
      else
        arg_values(no_args) = arg_values(no_args)//' '//flag%argument
      endif
    
    ! An argument which is preceded by '--', i.e. a long form option.
    elseif (flag%flag=='-') then
      no_args = no_args+1
      arg_keys(no_args) = flag%argument
      arg_values(no_args) = '<no_argument>'
    
    ! An argument which is preceded by '-', i.e. a flag.
    else
      ! Flags without arguments.
      do i=1,len(flags_without_arguments)
        if (flag%flag==slice(flags_without_arguments,i,i)) then
          no_args = no_args+1
          arg_keys(no_args) = long_flags_without_arguments(i)
          arg_values(no_args) = '<no_argument>'
          cycle do_args
        endif
      enddo
      
      ! Flags with arguments.
      do i=1,len(flags_with_arguments)
        if (flag%flag==slice(flags_with_arguments,i,i)) then
          no_args = no_args+1
          arg_keys(no_args) = long_flags_with_arguments(i)
          arg_values(no_args) = flag%argument
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
  arguments = make_dictionary(arg_keys(:no_args), arg_values(:no_args))
  
  ! Specify keywords for each mode.
  if (mode=='setup_harmonic') then
    keywords = setup_harmonic_keywords()
  elseif (mode/='') then
    call print_line('Error: unrecognised mode: '//mode)
    call print_line('Call caesar -h for help.')
    stop
  endif
  
  ! Handle help calls.
  temp_string = item(arguments, 'help')
  if (temp_string/='') then
    if (mode=='' .and. temp_string=='<no_argument>') then
      call help(mode)
    elseif (mode=='') then
      call help(temp_string)
    elseif (temp_string=='<no_argument>') then
      call help(mode, mode, keywords)
    else
      call help(temp_string, mode, keywords)
    endif
    stop
  endif
  
  ! Check that a mode has been specified.
  if (mode=='') then
    call print_line('Error: no mode specified. &
                    &Mode must be specified before any other arguments. &
                    &Call caesar -h for help.')
    stop
  endif
  
  ! Check if interactive mode requested.
  interactive = item(arguments, 'interactive')/=''
  
  ! Set working directory and input filename.
  wd = item(arguments, 'working_directory')
  if (wd=='<no_argument>') then
    call print_line('Error: no argument given with --working_directory.')
    stop
  endif
  
  input_filename = item(arguments, 'input_file')
  if (input_filename=='<no_argument>') then
    call print_line('Error: no argument given with --input_file.')
    stop
  endif
  
  ! Process input file.
  if (input_filename=='') then
    input_filename = '<not_set>'
  endif
  
  if (interactive) then
    call print_line('')
    call print_line('Input file is currently set to: '//input_filename)
    call print_line('Press <Enter> to accept current setting.')
    call print_line('Where is the input file?')
    temp_string = read_line_from_user()
    if (temp_string/='') then
      input_filename = temp_string
    endif
  endif
  
  if (input_filename/='<not_set>') then
    input_file = read_lines(input_filename)
    deallocate(arg_keys,arg_values,stat=ialloc); call err(ialloc)
    no_args = 0
    allocate( arg_keys(size(input_file)),   &
            & arg_values(size(input_file)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(input_file)
      line = split(input_file(i))
      ! Ignore empty lines.
      if (size(line)==0) then
        cycle
      
      ! Ignore comment lines (those beginning with a '!').
      elseif (slice(line(1),1,1)=='!') then
        cycle
      endif
      
      no_args = no_args+1
      arg_keys(no_args) = line(1)
      arg_values(no_args) = '<no_argument>'
      
      ! Read argument.
      do j=2,size(line)
        if (slice(line(j),1,1)=='!') then
          exit
        
        elseif (arg_values(no_args)=='<no_argument>') then
          arg_values(no_args) = line(j)
        
        else
          arg_values(no_args) = arg_values(no_args)//' '//line(j)
        endif
      enddo
    enddo
    
    ! Check for duplicate arguments.
    do i=1,no_args
      do j=i+1,no_args
        if (arg_keys(i)==arg_keys(j)) then
          call print_line('Error: argument '//arg_keys(i)// &
             & ' specified more than once in file '//input_filename//'.')
          stop
        endif
      enddo
    enddo
    
    ! Append file arguments to command line arguments.
    ! N.B. Command line arguments take priority over file arguments.
    arguments =  arguments &
               & // make_dictionary(arg_keys(:no_args), arg_values(:no_args))
  endif
  
  ! Process working directory.
  if (wd=='') then
    wd = '.'
  endif
  
  if (interactive) then
    call print_line('')
    call print_line('Working directory is currently set to: '//wd)
    call print_line('Press <Enter> to accept current setting.')
    call print_line('Where is the working directory?')
    temp_string = read_line_from_user()
    if (temp_string/='') then
      wd = temp_string
    endif
  endif
  
  wd = format_path(wd, cwd)
  
  ! Get interactive arguments.
  if (interactive) then
    no_args = 0
    deallocate(arg_keys, arg_values, stat=ialloc); call err(ialloc)
    allocate( arg_keys(size(keywords)),   &
            & arg_values(size(keywords)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(keywords)
      temp_string = item(arguments, keywords(i)%keyword)
      if (temp_string=='') then
        temp_string = '<not_set>'
        call print_line('')
      else
        call print_line('')
        call print_line(keywords(i)%keyword//' is currently set to: '//&
           & temp_string)
        call print_line('Press <Enter> to accept current setting.')
      endif
      
      do
        call print_line('Please give a value for the keyword: '// &
           & keywords(i)%keyword//'.')
        call print_line(keywords(i)%helptext)
        temp_string_2 = read_line_from_user()
        if (temp_string_2/='' .or. temp_string/='<not_set>') then
          exit
        endif
      enddo
      
      if (temp_string_2/='') then
        no_args = no_args+1
        arg_keys(no_args) = keywords(i)%keyword
        arg_values(no_args) = temp_string_2
      endif
    enddo
    
    ! Prepend interactive arguments to file and command line arguments.
    ! N.B. Interactive arguments take prioriry over all other arguments.
    arguments =  make_dictionary(arg_keys(:no_args), arg_values(:no_args)) &
            & // arguments
  endif
  
  ! Write settings to file.
  call write_dictionary_file(arguments,wd//'/settings.txt')
  
  ! Run main program, depending on mode.
  if (mode == 'test') then
    call system_call(str('mkdir this_is_another_dir'))
    call write_dictionary_file(arguments,wd//'/settings.txt')
  ! Wrappers for top-level Fortran.
  elseif (mode == 'setup_harmonic') then
    call setup_harmonic(wd, arguments)
  elseif (mode == 'run_harmonic') then
    call run_harmonic(wd,cwd)
  elseif (mode == 'lte_harmonic') then
    call lte_harmonic(wd)
  elseif (mode == 'setup_quadratic') then
    call setup_quadratic(wd,cwd)
  elseif (mode == 'run_quadratic') then
    call run_quadratic(wd,cwd)
  elseif (mode == 'anharmonic') then
    call anharmonic(wd)
  elseif (mode == 'bs_quadratic') then
    call bs_quadratic(wd)
  ! Wrappers for subsidiary Fortran 
  elseif (mode == 'calculate_gap') then
    call calculate_gap()
  elseif (mode == 'hartree_to_eV') then
    call hartree_to_eV()
  ! Wrappers for testing modules
  elseif (mode == 'test_copy_harmonic') then
    call test_copy_harmonic(wd,cwd)
  elseif (mode == 'test_lte') then
    call test_lte(wd,cwd)
  elseif (mode == 'test_copy_quadratic') then
    call test_copy_quadratic(wd,cwd)
  ! unrecognised mode.
  else
    call print_line('Error: Unrecognised mode: '//mode)
  endif
end program

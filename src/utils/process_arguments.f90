! ======================================================================
! Processes arguments from command line, file and interactive input.
! Stops with an error message if unexpected input is given.
! ======================================================================
module process_arguments_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! An interface for the main subroutines of Caesar, each of which takes a
  !    dictionary of arguments and returns nothing.
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
function process_arguments(args,keywords_in) result(arguments)
  use dictionary_module
  use help_module
  implicit none
  
  type(String),      intent(in) :: args(:)
  type(KeywordData), intent(in) :: keywords_in(:)
  type(Dictionary)              :: arguments
  
  type(KeywordData), allocatable :: universal_keywords(:)
  type(KeywordData), allocatable :: keywords(:)
  
  ! Flags.
  type(String)          :: flags_without_arguments
  type(String)          :: flags_with_arguments
  type(CommandLineFlag) :: flag
  
  ! Whether or not interactive mode is requested.
  logical :: interactive
  
  ! Temporary variables.
  integer      :: i,j,ialloc
  type(String) :: temp_string
  type(String) :: flags
  type(String) :: keyword
  logical      :: mode_found
  
  ! --------------------------------------------------
  ! Construct empty dictionary from keywords.
  ! --------------------------------------------------
  
  ! Concatenate input keywords with universal keywords ('help' etc.)
  universal_keywords = make_universal_keywords()
  allocate( keywords(size(keywords_in)+size(universal_keywords)), &
          & stat=ialloc); call err(ialloc)
  keywords(1:size(universal_keywords)) = universal_keywords
  keywords(size(universal_keywords)+1:) = keywords_in
  
  ! Check that keywords with default_keyword reference extant keywords.
  do_i : do i=1,size(keywords)
    if (keywords(i)%default_keyword/='') then
      do j=1,size(keywords)
        if (keywords(j)%keyword == keywords(i)%default_keyword) then
          cycle do_i
        endif
      enddo
      call print_line('Error: default keyword '//keywords(i)%default_keyword &
         & //'is not a keyword')
      call err()
    endif
  enddo do_i
  
  ! Parse allowed flags.
  flags_without_arguments = ''
  flags_with_arguments = ''
  do i=1,size(keywords)
    if (keywords(i)%flag/=' ') then
      if (keywords(i)%is_boolean) then
        flags_without_arguments = flags_without_arguments//keywords(i)%flag
      else
        flags_with_arguments = flags_with_arguments//keywords(i)%flag
      endif
    endif
  enddo
  
  ! Check that there are no duplicate flags.
  flags = flags_without_arguments//flags_with_arguments
  do i=1,len(flags)
    do j=i+1,len(flags)
      if (slice(flags,i,i)==slice(flags,j,j)) then
        call print_line('Error: the flag '//slice(flags,i,i)//' refers to two &
           &separate keywords.')
        call err()
      endif
    enddo
  enddo
  
  ! Convert keywords into Dictionary.
  arguments = keywords
  
  ! --------------------------------------------------
  ! Parse command line arguments.
  ! --------------------------------------------------
  keyword = ''
  mode_found = .false.
  do
    flag = get_flag(args, flags_without_arguments, flags_with_arguments)
    
    ! Check if the argument is preceded by '-' or '--'
    if (flag%flag==' ') then
      ! An argument which is neither a '-' flag or a '--' keyword.
      if (flag%argument==' ') then
        ! The end of the command line arguments
        exit
      elseif (keyword=='') then
        ! No '-' flag or '--' keyword has passed. This should be the mode.
        if (mode_found) then
          call print_line('Error: Caesar only takes one non-keyword &
             &argument. all other arguments should be preceded by "-" flags &
             &or "--" keywords.')
          stop
        else
          mode_found = .true.
        endif
      else
        ! Append this to the value of the active keyword.
        call print_line('arg: "'//flag%argument//'"')
        call arguments%append_value(keyword, flag%argument)
      endif
    elseif (flag%flag=='-') then
      ! An argument which is a '--' keyword.
      ! Get the keyword, and set its value to ''.
      keyword = lower_case(flag%argument)
    else
      ! An argument which is a '-' flag.
      ! Get the keyword corresponding to the flag, and set its value.
      keyword = arguments%flag_to_keyword(flag%flag)
      call arguments%set_value(keyword,flag%argument)
    endif
  enddo
  
  ! --------------------------------------------------
  ! Check if interactive mode is requested.
  ! --------------------------------------------------
  interactive = arguments%is_set('interactive')
  
  ! --------------------------------------------------
  ! Check if help is requested.
  ! --------------------------------------------------
  if (arguments%is_set('help')) then
    return
  endif
  
  if (interactive) then
    call print_line('')
    call print_line('Is help requested? Press <Enter> to skip or enter any &
       &value for help.')
    if (read_line_from_user()/='') then
      call print_line('Please enter a keyword for further information, or &
         &press <Enter> for information about accepted keywords.')
      call arguments%set_value('help',read_line_from_user())
      return
    endif
  endif
  
  ! --------------------------------------------------
  ! Process input file, if requested.
  ! --------------------------------------------------
  if (interactive) then
    if (arguments%is_set('input_file')) then
      call print_line('Current settings are to read arguments from file '// &
         & arguments%value('input_file'))
      call print_line('Please press <Enter> to confirm this, or enter any &
         &value to change this setting.')
      if (read_line_from_user()/='') then
        call arguments%unset('input_file')
      endif
    endif
    
    ! N.B. no elseif, because the above may unset input_file.
    if (.not. arguments%is_set('input_file')) then
      call print_line('Should further arguments be read from a file?')
      call print_line('Please enter a file path, or press <Enter> to skip.')
      temp_string = read_line_from_user()
      if (temp_string/='') then
        call arguments%set_value('input_file', temp_string)
      endif
    endif
  endif
  
  if (arguments%is_set('input_file')) then
    if (interactive) then
      call print_line('Please note that settings given on the command line &
         &override those read from file.')
    endif
    call arguments%read_file( arguments%value('input_file'), &
                       & only_set_if_not_set = .true. )
  endif
  
  ! --------------------------------------------------
  ! Get interactive arguments, if requested.
  ! --------------------------------------------------
  if (interactive) then
    do i=4,size(keywords) ! Skips 'interactive', 'help' and 'input_file'.
      keyword = keywords(i)%keyword
      call print_line('')
      call print_line(keywords(i)%helptext)
      if (arguments%is_set(keyword)) then
        if (keywords(i)%is_boolean) then
          ! Set boolean.
          call print_line(keyword//' is currently set.')
          call print_line('Please press <Enter> to leave it set, or enter any &
             &value to unset it.')
          if (read_line_from_user()/='') then
            call arguments%unset(keyword)
          endif
        else
          ! Set non-boolean.
          call print_line(keyword//' is currently set to '// &
             & arguments%value(keyword))
          call print_line('Please press <Enter> to accept this value, or &
             &enter the value to be set.')
          temp_string = read_line_from_user()
          if (temp_string/='') then
            call arguments%set_value(keyword, temp_string)
          endif
        endif
      else
        if (keywords(i)%is_boolean) then
          ! Unset boolean.
          call print_line(keyword//' is currently unset.')
          call print_line('Please press <Enter> to leave it unset, or enter &
             &any value to set it.')
          if (read_line_from_user()/='') then
            call arguments%set(keyword)
          endif
        elseif (keywords(i)%is_optional) then
          ! Unset optional argument.
          if (keywords(i)%default_value /= '') then
            call print_line(keyword//' defaults to the value '// &
               & keywords(i)%default_value)
            call print_line('Please give a value for this keyword, or press &
               &<Enter> to accept the default.')
          elseif (keywords(i)%default_keyword /= '') then
            call print_line(keyword//' defaults to the same value as &
               &keyword '//keywords(i)%default_keyword)
            call print_line('Please give a value for this keyword, or press &
               &<Enter> to accept the default.')
          else
            call print_line(keyword//' is optional.')
            call print_line('Please give a value for this keyword, or press &
               &<Enter> to leave it unset.')
          endif
          
          temp_string = read_line_from_user()
          if (temp_string/='') then
            call arguments%set_value(keyword, temp_string)
          endif
        else
          ! Unset non-optional argument. Loop until it is set by user.
          do
            call print_line('Please give a value for this keyword.')
            temp_string = read_line_from_user()
            if (temp_string/='') then
              call arguments%set_value(keyword, temp_string)
              exit
            endif
          enddo
        endif
      endif
    enddo
  endif
  
  ! --------------------------------------------------
  ! Set all defaults.
  ! Convert all paths to absolute format (from /).
  ! Check that all non-optional keywords have been set.
  ! --------------------------------------------------
  call arguments%process_and_check_inputs()
end function
end module

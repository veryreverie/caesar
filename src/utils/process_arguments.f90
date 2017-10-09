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
  use keyword_module
  implicit none
  
  type(String),      intent(in) :: args(:)
  type(KeywordData), intent(in) :: keywords_in(:)
  type(Dictionary)              :: arguments
  
  type(KeywordData), allocatable :: universal_keywords(:)
  type(KeywordData), allocatable :: keywords(:)
  
  ! Flags.
  type(String)          :: flags_without_arguments
  type(String)          :: flags_with_arguments
  type(CommandLineFlag) :: this
  
  ! Whether or not interactive mode is requested.
  logical :: interactive
  
  ! Temporary variables.
  integer      :: i,j,ialloc
  type(String) :: temp_string
  type(String) :: flags
  type(String) :: keyword
  logical      :: boolean_flag
  logical      :: mode_found
  type(String) :: default_keyword
  
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
    default_keyword = keywords(i)%defaults_to_keyword()
    if (default_keyword/='') then
      do j=1,size(keywords)
        if (keywords(j)%keyword == default_keyword) then
          cycle do_i
        endif
      enddo
      call print_line('Code Error: default keyword "'//default_keyword// &
         & '" is not a keyword')
      call err()
    endif
  enddo do_i
  
  ! Parse allowed flags.
  flags_without_arguments = ''
  flags_with_arguments = ''
  do i=1,size(keywords)
    if (keywords(i)%has_flag()) then
      if (keywords(i)%flag_takes_argument()) then
        flags_with_arguments = flags_with_arguments//keywords(i)%flag()
      else
        flags_without_arguments = flags_without_arguments//keywords(i)%flag()
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
    this = get_flag(args, flags_without_arguments, flags_with_arguments)
    
    ! ------------------------------
    ! Check if this is the first argument, which should be the mode.
    ! ------------------------------
    if (this%flag==' ' .and. this%argument/=' ' .and. keyword=='') then
      if (mode_found) then
        call print_line('Error: Caesar only takes one non-keyword &
           &argument. all other arguments should be preceded by "-" for flags &
           &or "--" for keywords.')
        stop
      endif
      mode_found = .true.
    endif
    
    ! ------------------------------
    ! Check the previous keyword has been given a value.
    ! ------------------------------
    ! Only required when this is a flag, keyword, or the end of args.
    if (keyword/='' .and. (this%flag/=' ' .or. this%argument==' ')) then
      if (.not. arguments%is_set(keyword)) then
        if (boolean_flag) then
          call arguments%set(keyword,'true')
        elseif (keyword=='help') then
          call arguments%set(keyword,'')
          return
        else
          call print_line('Error: keyword '//keyword//' has been given &
             &without a value on the command line.')
          stop
        endif
      endif
    endif
      
    ! ------------------------------
    ! Check if this is a keyword (preceeded by '--').
    ! ------------------------------
    if (this%flag=='-') then
      boolean_flag = .false.
      keyword = lower_case(this%argument)
    endif
    
    ! ------------------------------
    ! Check if this is a flag (preceeded by '-' or another flag).
    ! ------------------------------
    if (this%flag/='-' .and. this%flag/=' ') then
      if (index(char(flags_without_arguments),this%flag)/=0) then
        boolean_flag = .true.
      else
        boolean_flag = .false.
      endif
      keyword = arguments%flag_to_keyword(this%flag)
      if (this%argument/='') then
        call arguments%set(keyword,this%argument)
      endif
    endif
    
    ! ------------------------------
    ! Check if this is neither a flag nor a keyword.
    ! ------------------------------
    ! Append this to the value of the active keyword.
    if (this%flag==' ' .and. this%argument/=' ' .and. keyword/='') then
      if (arguments%is_set(keyword)) then
        call arguments%append_to_value(keyword, ' '//this%argument)
      else
        call arguments%set(keyword, this%argument)
      endif
    endif
    
    ! ------------------------------
    ! If none of the above are true, this should be the end of args.
    ! ------------------------------
    if (this%flag==' ' .and. this%argument==' ') then
      exit
    endif
  enddo
  
  ! --------------------------------------------------
  ! Check if interactive mode is requested.
  ! --------------------------------------------------
  if (arguments%is_set('interactive')) then
    interactive = lgcl(arguments%value('interactive'))
  else
    interactive = .false.
  endif
  
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
      call arguments%set('help',read_line_from_user())
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
        call arguments%set('input_file', temp_string)
      endif
    endif
  endif
  
  if (arguments%is_set('input_file')) then
    if (interactive) then
      call print_line('Please note that settings given on the command line &
         &override those read from file.')
    endif
    call arguments%read_file( arguments%value('input_file'), &
                            & only_update_if_unset = .true. )
  endif
  
  ! --------------------------------------------------
  ! Get interactive arguments, if requested.
  ! --------------------------------------------------
  if (interactive) then
    call arguments%set_interactively()
  endif
  
  ! --------------------------------------------------
  ! Set all defaults.
  ! Convert all paths to absolute format (from /).
  ! Check that all non-optional keywords have been set.
  ! --------------------------------------------------
  call arguments%process_and_check_inputs()
end function
end module

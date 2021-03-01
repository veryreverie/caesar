submodule (caesar_dictionary_module) caesar_dictionary_submodule
  use caesar_arguments_module
contains

module procedure new_Dictionary_KeywordDatas
  type(String), allocatable :: exclusives_i(:)
  type(String), allocatable :: exclusives_k(:)
  
  integer :: i,j,k,l
  
  this%keywords_ = [common_keywords(), keywords]
  
  ! Check that if keyword a is exclusive with keyword b then keyword b is also
  !    exclusive with keyword a.
  do i=1,size(this%keywords_)
    exclusives_i = this%keywords_(i)%exclusive_with()
    do j=1,size(exclusives_i)
      if (.not. any(this%keywords_%keyword()==exclusives_i(j))) then
        call print_line(CODE_ERROR//': keyword '// &
           & this%keywords_(i)%keyword()//' is mutually exclusive with &
           &keyword '//exclusives_i(j)//', which is not a keyword.')
        call err()
      endif
      
      k = this%index(exclusives_i(j))
      exclusives_k = this%keywords_(k)%exclusive_with()
      if (.not. any([( this%keywords_(i)%keyword()==exclusives_k(l), &
                     & l=1,                                          &
                     & size(exclusives_k)                            )])) then
        call print_line(CODE_ERROR//': mutually exclusive keyword lists &
           &inconsistent between keywords '//this%keywords_(i)%keyword()//' &
           &and '//this%keywords_(k)%keyword()//'.')
        call err()
      endif
    enddo
  enddo
end procedure

module procedure new_Dictionary_arguments
  type(KeywordData), allocatable :: keywords(:)
  
  ! Flags.
  type(String)          :: flags_without_arguments
  type(String)          :: flags_with_arguments
  type(CommandLineFlag) :: this
  
  ! Whether or not interactive mode is requested.
  logical :: interactive
  
  ! Temporary variables.
  integer      :: i,j
  type(String) :: temp_string
  type(String) :: flags
  type(String) :: keyword
  logical      :: boolean_flag
  logical      :: mode_found
  type(String) :: default_keyword
  
  ! --------------------------------------------------
  ! Construct empty dictionary from keywords.
  ! --------------------------------------------------
  
  ! Concatenate input keywords with common keywords ('help' etc.)
  keywords = [common_keywords(), keywords_in]
  
  ! Check that keywords with default_keyword reference extant keywords.
  do_i : do i=1,size(keywords)
    default_keyword = keywords(i)%defaults_to_keyword()
    if (default_keyword/='') then
      do j=1,size(keywords)
        if (keywords(j)%keyword() == default_keyword) then
          cycle do_i
        endif
      enddo
      call print_line(CODE_ERROR//': default keyword "'//default_keyword// &
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
        call print_line(CODE_ERROR//': the flag '//slice(flags,i,i)//' refers &
           &to two separate keywords.')
        call err()
      endif
    enddo
  enddo
  
  ! Convert keywords into Dictionary.
  ! N.B. the Dictionary() constructor will automatically prepend common
  !    keywords to keywords_in.
  arguments = Dictionary(keywords_in)
  
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
        call print_line(ERROR//': Caesar only takes one non-keyword &
           &argument. all other arguments should be preceded by "-" for flags &
           &or "--" for keywords.')
        call quit()
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
          call print_line(ERROR//': keyword '//keyword//' has been given &
             &without a value on the command line.')
          call quit()
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
end procedure

module procedure size_Dictionary
  output = size(this%keywords_)
end procedure

module procedure index_Dictionary_String
  output = first(this%keywords_%keyword() == keyword, default=0)
  
  if (output==0) then
    call print_line(ERROR//': unexpected keyword: '//keyword//'.')
    call quit()
  endif
end procedure

module procedure index_Dictionary_character
  output = this%index(str(keyword))
end procedure

module procedure index_by_flag_Dictionary_character
  output = first(this%keywords_,flag_matches,default=0)
  
  if (output==0) then
    call print_line(ERROR//': unexpected flag: '//flag//'.')
    call quit()
  endif
contains
  ! Lambda for checking if a flag matches.
  ! Captures flag from the module function.
  function flag_matches(input) result(output) 
    class(*), intent(in) :: input
    logical              :: output
    select type(input); class is(KeywordData)
      if (input%has_flag()) then
        output = input%flag()==flag
      else
        output = .false.
      endif
    end select
  end function
end procedure

module procedure index_by_flag_Dictionary_String
  output = this%index_by_flag(char(slice(flag,1,1)))
end procedure

module procedure flag_to_keyword_Dictionary_character
  output = this%keywords_(this%index_by_flag(flag))%keyword()
end procedure

module procedure flag_to_keyword_Dictionary_String
  output = this%flag_to_keyword(char(flag))
end procedure

module procedure is_set_Dictionary_character
  output = this%keywords_(this%index(keyword))%is_set()
end procedure

module procedure is_set_Dictionary_String
  output = this%is_set(char(keyword))
end procedure

module procedure value_Dictionary_character
  output = this%keywords_(this%index(keyword))%value()
end procedure

module procedure value_Dictionary_String
  output = this%value(char(keyword))
end procedure

module procedure python_arguments
  integer :: i,ialloc
  
  allocate(output(0), stat=ialloc); call err(ialloc)
  do i=1,size(this%keywords_)
    if ( this%keywords_(i)%is_set() .and.   &
       & this%keywords_(i)%pass_to_python() ) then
      output = [ output,                                                     &
               & this%keywords_(i)%keyword()//' '//this%keywords_(i)%value() ]
    endif
  enddo
end procedure

module procedure unset_Dictionary_character
  call this%keywords_(this%index(keyword))%unset()
end procedure

module procedure unset_Dictionary_String
  call this%unset(char(keyword))
end procedure

module procedure set_Dictionary_character_character
  call this%keywords_(this%index(keyword))%set(value, only_update_if_unset) 
end procedure

module procedure set_Dictionary_character_String
  call this%set(keyword, char(value), only_update_if_unset)
end procedure

module procedure set_Dictionary_String_character
  call this%set(char(keyword), value, only_update_if_unset)
end procedure

module procedure set_Dictionary_String_String
  call this%set(char(keyword), char(value), only_update_if_unset)
end procedure

module procedure set_Dictionary_Dictionary
  type(String) :: keyword
  
  integer :: i
  
  do i=1,size(this%keywords_)
    keyword = this%keywords_(i)%keyword()
    if (keyword/='interactive' .and. keyword/='input_file') then
      if (any(keyword==that%keywords_%keyword())) then
        if (that%is_set(keyword)) then
          call this%set(keyword, that%value(keyword))
        endif
      endif
    endif
  enddo
end procedure

module procedure append_to_value_Dictionary_character_character
  call this%keywords_(this%index(keyword))%append(value)
end procedure

module procedure append_to_value_Dictionary_character_String
  call this%append_to_value(keyword, char(value))
end procedure

module procedure append_to_value_Dictionary_String_character
  call this%append_to_value(char(keyword), value)
end procedure

module procedure append_to_value_Dictionary_String_String
  call this%append_to_value(char(keyword), char(value))
end procedure

module procedure write_file_Dictionary_character
  type(OFile) :: dictionary_file
  
  logical :: to_write
  integer :: max_length
  
  integer :: i
  
  ! Check that there are any settings to write.
  ! If not, return without creating a blank file.
  ! Also get the length of the longest keyword, for formatting purposes.
  to_write = .false.
  max_length = 0
  do i=1,size(this)
    if (.not. this%keywords_(i)%allowed_in_file()) then
      cycle
    elseif (.not. this%keywords_(i)%is_set()) then
      cycle
    elseif (len(this%keywords_(i)%value())==0) then
      cycle
    else
      to_write = .true.
      max_length = max(max_length, len(this%keywords_(i)%keyword()))
    endif
  enddo
  
  if (to_write) then
    dictionary_file = OFile(filename)
    do i=1,size(this)
      if (.not. this%keywords_(i)%allowed_in_file()) then
        cycle
      elseif (.not. this%keywords_(i)%is_set()) then
        cycle
      elseif (len(this%keywords_(i)%value())==0) then
        cycle
      else
        call dictionary_file%print_line(                              &
           & this%keywords_(i)%keyword() //                           &
           & spaces(max_length-len(this%keywords_(i)%keyword())+1) // &
           & this%keywords_(i)%value())
      endif
    enddo
  endif
end procedure

module procedure write_file_Dictionary_String
  call this%write_file(char(filename))
end procedure

module procedure read_file_Dictionary_character
  ! Files.
  type(IFile) :: dictionary_file
  
  ! The directory where the file is.
  type(String) :: file_directory
  
  ! Temporary variables.
  integer                   :: i,j
  type(String), allocatable :: line(:)
  logical                   :: only_if_unset ! = only_update_if_unset
  
  ! Set default value of only_update_if_unset.
  if (present(only_update_if_unset)) then
    only_if_unset = only_update_if_unset
  else
    only_if_unset = .false.
  endif
  
  ! Parse the directory where the file is.
  line = split_line(format_path(filename), '/')
  if (size(line)<2) then
    call err()
  else
    file_directory = '/'//join(line(1:size(line)-1),'/')
  endif
  
  ! Read file.
  dictionary_file = IFile(filename)
  
  ! Process file.
  ! Each line is expected to be of the form '  key  value  ! comments  '
  do i=1,size(dictionary_file)
    ! Strip leading and trailing spaces.
    line = [trim(dictionary_file%line(i))]
    
    ! Ignore blank lines and comment lines (those beginning with a '!').
    if (len(line(1))==0) then
      cycle
    elseif (slice(line(1),1,1)=='!') then
      cycle
    endif
    
    ! Strip out anything after a comment character (!).
    line = split_line(line(1), '!')
    line = split_line(line(1))
    
    ! Find keyword in arguments.
    j = this%index(lower_case(line(1)))
    if (.not. this%keywords_(j)%allowed_in_file()) then
      call print_line(ERROR//': the keyword '//this%keywords_(j)%keyword()// &
         & ' should not appear in input files.')
      call err()
    elseif (only_if_unset .and. this%keywords_(j)%is_set()) then
      call print_line(WARNING//': the keyword '//    &
         & this%keywords_(j)%keyword()         //    &
         & ' has been specified in multiple places.' )
    else
      if (size(line)==1) then
        call print_line(ERROR//': the keyword '//  &
           & this%keywords_(j)%keyword()       //  &
           & 'has been specified without a value.' )
        call quit()
      else
        call this%keywords_(j)%set( value                = join(line(2:)), &
                                  & only_update_if_unset = only_if_unset,  &
                                  & working_directory    = file_directory  )
      endif
    endif
  enddo
end procedure

module procedure read_file_Dictionary_String
  call this%read_file(char(filename), only_update_if_unset)
end procedure

module procedure set_interactively_Dictionary
  type(String), allocatable :: exclusives(:)
  
  integer :: i
  
  do i=1,size(this)
    ! Check if this keyword can be set interactively.
    if (this%keywords_(i)%can_be_interactive()) then
      ! Check if this keyword is mutually exclusive with a keyword which has
      !    already been set.
      exclusives = this%keywords_(i)%exclusive_with()
      exclusives = exclusives(filter(this%index(exclusives)<i))
      if (.not.any([(this%is_set(exclusives(i)),i=1,size(exclusives))])) then
        ! This keyword is not mutually exclusive with a keyword which has
        !    already been set; set it interactively.
        call this%keywords_(i)%set_interactively()
      else
        ! This keyword is mutually exclusive with a keyword which has
        !    already been set; unset it.
        call this%keywords_(i)%unset()
      endif
    endif
  enddo
end procedure

module procedure process_and_check_inputs_Dictionary
  type(String), allocatable :: exclusives(:)
  
  type(String) :: default_keyword
  
  integer :: i,j,k
  
  do_i : do i=1,size(this)
    if (this%keywords_(i)%is_set()) then
      exclusives = this%keywords_(i)%exclusive_with()
      j = first([(this%is_set(exclusives(k)),k=1,size(exclusives))], default=0)
      if (j/=0) then
        call print_line(ERROR//': the keywords '//this%keywords_(i)%keyword() &
           &//' and '//exclusives(j)//' are mutually exclusive but have both &
           &been set.')
        call quit()
      endif
    else
      default_keyword = this%keywords_(i)%defaults_to_keyword()
      if (default_keyword == '') then
        call this%keywords_(i)%set_default()
      else
        do j=1,size(this)
          if (this%keywords_(j)%keyword()==default_keyword) then
            if (this%keywords_(j)%is_set()) then
              call this%keywords_(i)%set(this%keywords_(j)%value())
            endif
            cycle do_i
          endif
        enddo
        call err()
      endif
    endif
  enddo do_i
  
  do i=1,size(this)
    call this%keywords_(i)%process_and_check()
  enddo
end procedure

module procedure call_caesar_String_Dictionary
  type(String) :: command_line_arguments
  
  integer :: i
  
  command_line_arguments = mode
  do i=1,size(arguments%keywords_)
    if (arguments%keywords_(i)%is_set()) then
      if (arguments%keywords_(i)%value()/='') then
        command_line_arguments = command_line_arguments                 // &
                               & ' '                                    // &
                               & '--'//arguments%keywords_(i)%keyword() // &
                               & ' '                                    // &
                               & arguments%keywords_(i)%value()
      endif
    endif
  enddo
  
  call call_caesar(command_line_arguments)
end procedure

module procedure call_caesar_character_Dictionary
  call call_caesar(str(mode), arguments)
end procedure
end submodule

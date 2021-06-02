submodule (caesar_keyword_module) caesar_keyword_submodule
  use caesar_arguments_module
contains

module procedure has_flag_KeywordData
  output = this%flag_type_ > 0
end procedure

module procedure flag_takes_argument_KeywordData
  output = this%flag_type_ == 2
end procedure

module procedure flag_KeywordData
  if (this%flag_type_==0) then
    call err()
  endif
  output = this%flag_
end procedure

module procedure set_flag_KeywordData
  if (flag_takes_arguments) then
    this%flag_type_ = 2
  else
    this%flag_type_ = 1
  endif
  this%flag_ = flag
end procedure

module procedure defaults_to_keyword_KeywordData
  if (this%default_type_/=3) then
    output = ''
  else
    output = this%default_
  endif
end procedure

module procedure set_default_KeywordData
  if (this%default_type_==2 .and. .not. this%is_set()) then
    call this%set(this%default_)
  endif
end procedure

module procedure unset_KeywordData
  this%is_set_ = .false.
end procedure

module procedure set_KeywordData_character
  logical :: only_update_if
  
  integer :: ialloc
  
  if (present(only_update_if_unset)) then
    only_update_if = only_update_if_unset
  else
    only_update_if = .false.
  endif
  
  if (.not. (only_update_if .and. this%is_set())) then
    this%is_set_ = .true.
    this%value_ = value
    
    if (present(working_directory)) then
      this%working_directory_ = working_directory
    else
      if (allocated(this%working_directory_)) then
        deallocate(this%working_directory_, stat=ialloc); call err(ialloc)
      endif
    endif
  endif
end procedure

module procedure set_KeywordData_String
  call this%set(char(value),only_update_if_unset,working_directory)
end procedure

module procedure append_KeywordData_character
  if (.not. this%is_set()) then
    call err()
  endif
  this%value_ = this%value_ // value
end procedure

module procedure append_KeywordData_String
  call this%append(char(value))
end procedure

module procedure keyword_KeywordData
  output = this%keyword_
end procedure

module procedure is_set_KeywordData
  output = this%is_set_
end procedure

module procedure value_KeywordData
  if (.not. this%is_set()) then
    call print_line(CODE_ERROR//': Requesting the value of keyword '// &
       & this%keyword()//', which has not been set.')
    call err()
  endif
  output = this%value_
end procedure

module procedure is_path_KeywordData
  output = this%is_path_
end procedure

module procedure allowed_in_file_KeywordData
  output = this%allowed_in_file_
end procedure

module procedure can_be_interactive_KeywordData
  output = this%can_be_interactive_
end procedure

module procedure pass_to_python_KeywordData
  output = this%pass_to_python_
end procedure

module procedure exclusive_with_KeywordData
  output = this%exclusive_with_
end procedure

module procedure new_KeywordData
  integer :: ialloc
  
  ! Check for incompatible arguments.
  if (present(is_optional)) then
    if ( present(default_value)   .or. &
       & present(default_keyword) .or. &
       & present(exclusive_with)       ) then
      call print_line(CODE_ERROR//': the argument "is_optional" should not be &
         &given for keywords which have defaults or which are mutually &
         &exclusive with other keywords.')
      call err()
    endif
  endif
  
  if (present(default_value) .and. present(default_keyword)) then
    call print_line(CODE_ERROR//': a keyword may not have a default value and &
       &default to another keyword.')
    call err()
  endif
  
  if (present(flag_without_arguments) .and. present(flag_with_arguments)) then
    call print_line(CODE_ERROR//': a keyword may not have two flags.')
    call err()
  endif
  
  ! Set properties.
  this%keyword_ = lower_case(keyword)
  this%helptext_ = helptext
  
  if (present(is_path)) then
    this%is_path_ = is_path
  else
    this%is_path_ = .false.
  endif
  
  if (present(allowed_in_file)) then
    this%allowed_in_file_ = allowed_in_file
  else
    this%allowed_in_file_ = .true.
  endif
  
  if (present(can_be_interactive)) then
    this%can_be_interactive_ = can_be_interactive
  else
    this%can_be_interactive_ = .true.
  endif
  
  if (present(flag_without_arguments)) then
    this%flag_type_ = 1
    this%flag_ = flag_without_arguments
  elseif (present(flag_with_arguments)) then
    this%flag_type_ = 2
    this%flag_ = flag_with_arguments
  else
    this%flag_type_ = 0
  endif
  
  if (present(pass_to_python)) then
    this%pass_to_python_ = pass_to_python
  else
    this%pass_to_python_ = .false.
  endif
  
  if (present(exclusive_with)) then
    if (size(exclusive_with)==0) then
      call print_line(ERROR//': exclusive_with has been specified, but no &
         &keywords given.')
      call err()
    endif
    this%exclusive_with_ = exclusive_with
  else
    allocate(this%exclusive_with_(0), stat=ialloc); call err(ialloc)
  endif
  
  ! Set default behaviour.
  if (present(is_optional)) then
    if (is_optional) then
      this%default_type_ = 1
    else
      this%default_type_ = 0
    endif
  elseif (present(default_value)) then
    this%default_type_ = 2
    this%default_ = default_value
  elseif (present(default_keyword)) then
    this%default_type_ = 3
    this%default_ = default_keyword
  elseif (present(exclusive_with)) then
    this%default_type_ = 1
  else
    this%default_type_ = 0
  endif
  
  ! Unset value.
  call this%unset()
end procedure

module procedure set_interactively
  type(String), allocatable :: exclusives(:)
  
  type(String) :: input
  
  call print_line('')
  call print_line(this%helptext_)
  exclusives = this%exclusive_with_
  if (size(exclusives)==1) then
    call print_line('This keyword is mutually exclusive with keyword '// &
       & exclusives(1)                                                // &
       & '. Only one of these keywords should be specified.'             )
  elseif (size(exclusives)>1) then
    call print_line('This keyword is mutually exclusive with keywords '// &
       & join(exclusives(:size(exclusives)-1),delimiter=', ')          // &
       & ' and '                                                       // &
       & exclusives(size(exclusives))                                  // &
       & '. Only one of these keywords should be specified.'              )
  endif
  if (this%is_set()) then
    call print_line(this%keyword_//' currently has the value "'//this%value() &
       & //'".')
    call print_line('Please press <Enter> to accept this value, or enter &
       &anything for other options.')
    input = read_line_from_user()
    if (input/='') then
      call this%unset()
      call this%set_interactively()
    endif
  else
    if (this%default_type_==0) then
      call print_line(this%keyword_//' is unset and has no default.')
      do while (.not. this%is_set())
        call print_line('Please enter a value.')
        input = read_line_from_user()
        if (input/='') then
          call this%set(input)
        endif
      enddo
    else
      if (this%default_type_==1) then
        call print_line(this%keyword_//' is unset and is optional.')
        call print_line('Please press <Enter> to leave it unset, or enter a &
           &value.')
      elseif (this%default_type_==2) then
        call print_line(this%keyword_//' is unset, and will default to the &
           &value "'//this%default_//'".')
        call print_line('Please press <Enter> to accept this value, or enter &
           &a value.')
      elseif (this%default_type_==3) then
        call print_line(this%keyword_//' is unset, and will default to the &
           &keyword "'//this%default_//'".')
        call print_line('Please press <Enter> to accept this default, or &
           &enter a value.')
      endif
      
      input = read_line_from_user()
      if (input/='') then
        call this%set(input)
      endif
    endif
  endif
end procedure

module procedure process_and_check
  if (this%default_type_==0 .and. .not. this%is_set()) then
    call print_line(ERROR//': the keyword '//this%keyword_//' has not been &
       &set. This keyword is not optional.')
    call quit()
  endif
  
  if (this%is_path_ .and. this%is_set()) then
    this%value_ = format_path(this%value_, this%working_directory_)
  endif
end procedure

module procedure print_help
  type(String), allocatable :: helptext(:)
  type(String), allocatable :: exclusives(:)
  integer                   :: i
  
  ! Find the first instance of the keyword in the helptext,
  !    and colour it white.
  helptext = split_line(this%helptext_)
  i = first(helptext==this%keyword_, default=0)
  if (i/=0) then
    helptext(i) = colour(helptext(i), 'white')
  endif
  
  call print_line('')
  call print_line(join(helptext))
  
  exclusives = [( colour(this%exclusive_with_(i),'white'), &
                & i=1,                                     &
                & size(this%exclusive_with_)               )]
  if (size(exclusives)==1) then
    call print_line('This keyword is mutually exclusive with keyword '// &
       & exclusives(1)                                                // &
       & '. Only one of these keywords should be specified.'             )
  elseif (size(exclusives)>1) then
    call print_line('This keyword is mutually exclusive with keywords '// &
       & join(exclusives(:size(exclusives)-1),delimiter=', ')          // &
       & ' and '                                                       // &
       & exclusives(size(exclusives))                                  // &
       & '. Only one of these keywords should be specified.'              )
  endif
  
  if (this%default_type_==0) then
    call print_line('This keyword is non-optional.')
  elseif (this%default_type_==1) then
    call print_line('This keyword is optional but has no default.')
  elseif (this%default_type_==2) then
    call print_line('This keyword has a default value of "'//this%default_// &
       & '".')
  elseif (this%default_type_==3) then
    call print_line('This keyword defaults to the same value as keyword "'// &
       & this%default_//'".')
  else
    call err()
  endif
end procedure
end submodule

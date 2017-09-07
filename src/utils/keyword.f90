! ======================================================================
! Keywords for input arguments.
! ======================================================================
module keyword_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  ! A keyword and its helptext.
  public :: KeywordData
  
  ! Bundles a keyword and its helptext into a KeywordData.
  public :: make_keyword
  
  ! Keywords used by all routines.
  public :: make_universal_keywords
  
  ! --------------------------------------------------
  ! A keyword which accepts input values.
  ! --------------------------------------------------
  type KeywordData
    ! Keyword and flag (alternative one-character keyword).
    type(String) :: keyword
    character(1) :: flag
    
    ! Helptext for help calls.
    type(String) :: helptext
    
    ! Defaults.
    type(String) :: default_value
    type(String) :: default_keyword
    
    ! Properties of the value.
    logical      :: is_boolean
    logical      :: is_optional
    logical      :: is_path
    logical      :: allowed_in_file
    
    ! The value itself.
    logical      :: is_set
    logical      :: is_set_with_value
    type(String) :: value
  end type
contains

! ----------------------------------------------------------------------
! Takes a keyword, helptext and options and returns a KeywordData.
! ----------------------------------------------------------------------
! default_value is the value which will be defaulted to if the keyword is not
!    specified. Defaults to not set.
! default_keyword is the keyword whose value will be copied if the keyword is
!    not specified. Defaults to not set.
! At most one of default_value and default_keyword should be specified.
! If is_boolean is .true., the keyword will not take an argument.
!    Boolean keywords should not have defaults specified.
!    is_boolean defaults to .false..
! is_optional specifies whether or not the keyword is required to be given.
!    is_optional defaults to .false. unless a default is given or is_boolean
!    is .true., where it defaults to and is required to be .true..
! is_path specifies that the value will be a path to a file or directory.
!    Paths are converted to absolute paths (from /) automatically.
!    Should not be specified if is_boolean is .true..
!    Defaults to .false..
! allowed_in_file specifies that the value may be read from or printed to file.
!    Defaults to .true..
! flag specifies the flag by which the keyword can be alternately called.
!    Defaults to ' '.
function make_keyword(keyword,helptext,default_value,default_keyword, &
   & is_boolean,is_optional,is_path,allowed_in_file,flag) result(this)
  implicit none
  
  character(*), intent(in)           :: keyword
  character(*), intent(in)           :: helptext
  character(*), intent(in), optional :: default_value
  character(*), intent(in), optional :: default_keyword
  logical,      intent(in), optional :: is_boolean
  logical,      intent(in), optional :: is_optional
  logical,      intent(in), optional :: is_path
  logical,      intent(in), optional :: allowed_in_file
  character(1), intent(in), optional :: flag
  type(KeywordData)                  :: this
  
  ! Set keyword; flag, if present; and helptext.
  this%keyword = lower_case(keyword)
  
  if (present(flag)) then
    this%flag = flag
  else
    this%flag = ' '
  endif
  this%helptext = helptext
  
  ! Set type of keyword.
  if (present(is_boolean)) then
    this%is_boolean = is_boolean
  else
    this%is_boolean = .false.
  endif
  
  if (present(is_optional)) then
    this%is_optional = is_optional
  else
    if (this%is_boolean) then
      this%is_optional = .true.
    elseif (present(default_value) .or. present(default_keyword)) then
      this%is_optional = .true.
    else
      this%is_optional = .false.
    endif
  endif
  
  if (present(is_path)) then
    this%is_path = is_path
  else
    this%is_path = .false.
  endif
  
  if (present(allowed_in_file)) then
    this%allowed_in_file = allowed_in_file
  else
    this%allowed_in_file = .true.
  endif
  
  if (this%is_boolean) then
    ! Check for conflicts with the keyword being boolean.
    if (this%is_path) then
      call print_line('Error: the keyword '//keyword//' cannot be boolean and &
         &be a path.')
      call err()
    elseif (.not. this%is_optional) then
      call print_line('Error: the keyword '//keyword//' cannot be boolean and &
         &not be optional.')
      call err()
    elseif (present(default_value)) then
      call print_line('Error: the keyword '//keyword//' cannot be boolean and &
         &have a default value.')
      call err()
    endif
    
    ! Set default value, if present.
    this%default_value = ''
    if (present(default_keyword)) then
      this%default_keyword = lower_case(default_keyword)
    else
      this%default_keyword = ''
    endif
  else
    ! Check for conflicts between defaults.
    if (present(default_value) .and. present(default_keyword)) then
      call print_line('Error: the keyword '//keyword//' cannot have both a &
         &default value and a default keyword set.')
      call err()
    elseif ( (present(default_value) .or. present(default_keyword)) .and. &
           & .not. this%is_optional) then
      call print_line('Error: the keyword '//keyword//' cannot have a default &
         &and not be optional.')
        call err()
    endif
    
    ! Set defaults, if present.
    if (present(default_value)) then
      this%default_value = default_value
      this%default_keyword = ''
    elseif (present(default_keyword)) then
      this%default_value = ''
      this%default_keyword = lower_case(default_keyword)
    else
      this%default_value = ''
      this%default_keyword = ''
    endif
  endif
  
  ! Set value to unset state.
  this%is_set = .false.
  this%is_set_with_value = .false.
end function

! ----------------------------------------------------------------------
! Universal keywords.
! ----------------------------------------------------------------------
function make_universal_keywords() result(keywords)
  implicit none
  
  type(KeywordData) :: keywords(4)
  
  keywords = [                                                                &
  & make_keyword( 'interactive',                                              &
  &               'interactive specifies whether or not keywords can be &
  &specified interactively.',                                                 &
  &               is_boolean=.true.,                                          &
  &               allowed_in_file=.false.,                                    &
  &               flag='i'),                                                  &
  & make_keyword( 'help',                                                     &
  &               'help requests helptext rather than running calculation.', &
  &               is_optional=.true.,                                         &
  &               allowed_in_file=.false.,                                    &
  &               flag='h'),                                                  &
  & make_keyword( 'input_file',                                               &
  &               'input_file specifies a file from which further settings &
  &will be read.',                                                            &
  &               is_optional=.true.,                                         &
  &               allowed_in_file=.false.,                                    &
  &               flag='f'),                                                  &
  & make_keyword( 'working_directory',                                        &
  &               'working_directory specifies the directory where all files &
  &and subsequent directories will be made.',                                 &
  &               default_value='.',                                          &
  &               allowed_in_file=.false.,                                    &
  &               flag='d') ]
end function
end module

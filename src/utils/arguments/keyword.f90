! ======================================================================
! Keywords for input arguments.
! ======================================================================
module keyword_submodule
  use precision_module
  use io_module
  implicit none
  
  private
  
  public :: KeywordData
  
  ! A keyword for getting input arguments from the user.
  type KeywordData
    ! Keyword.
    type(String) :: keyword
    
    ! Flag (alternative one-character keyword).
    ! 0=no flag,1=flag without arguments,2=flag which takes arguments.
    integer,      private :: flag_type_
    character(1), private :: flag_
    
    ! Helptext for help calls.
    type(String) :: helptext
    
    ! Properties of the value.
    logical :: is_path
    logical :: allowed_in_file
    logical :: can_be_interactive
    
    ! Defaults.
    ! 0=must be set,1=no default,2=default to value,3=default to keyword.
    integer,      private :: default_type_
    type(String), private :: default_
    
    ! The value itself.
    logical,      private :: is_set_
    type(String), private :: value_
  contains
    ! Flag-related procedures.
    procedure, public :: has_flag            => has_flag_KeywordData
    procedure, public :: flag_takes_argument => flag_takes_argument_KeywordData
    procedure, public :: flag                => flag_KeywordData
    procedure, public :: set_flag            => set_flag_KeywordData
    
    ! Default-related procedures.
    procedure, public :: defaults_to_keyword
    procedure, public :: set_default
    
    ! Setters.
    procedure, public  :: unset   => unset_KeywordData
    generic,   public  :: set     => set_KeywordData_character, &
                                   & set_KeywordData_String
    procedure, private ::            set_KeywordData_character
    procedure, private ::            set_KeywordData_String
    generic,   public  :: append  => append_KeywordData_character, &
                                   & append_KeywordData_String
    procedure, private ::            append_KeywordData_character
    procedure, private ::            append_KeywordData_String
    
    ! Getters.
    procedure, public :: is_set => is_set_KeywordData
    procedure, public :: value  => value_KeywordData
    
    ! Set interactively from user.
    procedure, public :: set_interactively
    
    ! Processes and checks value.
    procedure, public :: process_and_check
    
    ! Prints helptext.
    procedure, public :: print_help
  end type
  
  interface KeywordData
    module procedure new_KeywordData
  end interface
contains

! ----------------------------------------------------------------------
! Flag-related procedures.
! ----------------------------------------------------------------------
function has_flag_KeywordData(this) result(output)
  implicit none
  
  class(KeywordData), intent(in) :: this
  logical                        :: output
  
  output = this%flag_type_ > 0
end function

function flag_takes_argument_KeywordData(this) result(output)
  implicit none
  
  class(KeywordData), intent(in) :: this
  logical                        :: output
  
  output = this%flag_type_ == 2
end function

function flag_KeywordData(this) result(output)
  implicit none
  
  class(KeywordData), intent(in) :: this
  character(1)                   :: output
  
  if (this%flag_type_==0) then
    call err()
  endif
  output = this%flag_
end function

subroutine set_flag_KeywordData(this,flag,flag_takes_arguments)
  implicit none
  
  class(KeywordData), intent(inout) :: this
  character(1),       intent(in)    :: flag
  logical,            intent(in)    :: flag_takes_arguments
  
  if (flag_takes_arguments) then
    this%flag_type_ = 2
  else
    this%flag_type_ = 1
  endif
  this%flag_ = flag
end subroutine

! ----------------------------------------------------------------------
! Default-related procedures.
! ----------------------------------------------------------------------

! If the keyword defaults to another keyword, returns that keyword.
! Returns '' if the keyword does not default to a keyword.
function defaults_to_keyword(this) result(output)
  implicit none
  
  class(KeywordData), intent(in) :: this
  type(String)                   :: output
  
  if (this%default_type_/=3) then
    output = ''
  else
    output = this%default_
  endif
end function

! Sets an unset keyword to its default value.
! Does nothing if the keyword is set or has no default.
! Throws an error if the keyword defaults to another keyword.
subroutine set_default(this)
  implicit none
  
  class(KeywordData), intent(inout) :: this
  
  if (this%default_type_==2 .and. .not. this%is_set()) then
    call this%set(this%default_)
  endif
end subroutine

! ----------------------------------------------------------------------
! Setters.
! ----------------------------------------------------------------------
subroutine unset_KeywordData(this)
  implicit none
  
  class(KeywordData), intent(inout) :: this
  
  this%is_set_ = .false.
end subroutine

subroutine set_KeywordData_character(this,value,only_update_if_unset)
  implicit none
  
  class(KeywordData), intent(inout)        :: this
  character(*),       intent(in)           :: value
  logical,            intent(in), optional :: only_update_if_unset
  
  logical :: only_update_if
  
  if (present(only_update_if_unset)) then
    only_update_if = only_update_if_unset
  else
    only_update_if = .false.
  endif
  
  if (.not. (only_update_if .and. this%is_set())) then
    this%is_set_ = .true.
    this%value_ = value
  endif
end subroutine

subroutine set_KeywordData_String(this,value,only_update_if_unset)
  implicit none
  
  class(KeywordData), intent(inout)        :: this
  type(String),       intent(in)           :: value
  logical,            intent(in), optional :: only_update_if_unset
  
  if (present(only_update_if_unset)) then
    call this%set(char(value),only_update_if_unset)
  else
    call this%set(char(value))
  endif
end subroutine

subroutine append_KeywordData_character(this,value)
  implicit none
  
  class(KeywordData), intent(inout)        :: this
  character(*),       intent(in)           :: value
  
  if (.not. this%is_set()) then
    call err()
  endif
  this%value_ = this%value_ // value
end subroutine

subroutine append_KeywordData_String(this,value)
  implicit none
  
  class(KeywordData), intent(inout)        :: this
  type(String),       intent(in)           :: value
  
  call this%append(char(value))
end subroutine

! ----------------------------------------------------------------------
! Getters.
! ----------------------------------------------------------------------
function is_set_KeywordData(this) result(output)
  implicit none
  
  class(KeywordData), intent(in) :: this
  logical                        :: output
  
  output = this%is_set_
end function

function value_KeywordData(this) result(output)
  implicit none
  
  class(KeywordData), intent(in) :: this
  type(String)                   :: output
  
  if (.not. this%is_set()) then
    call err()
  endif
  output = this%value_
end function

! ----------------------------------------------------------------------
! Takes a keyword, helptext and options and returns a KeywordData.
! ----------------------------------------------------------------------
! default_value is the value which will be defaulted to if the keyword is not
!    specified. Defaults to not set.
! default_keyword is the keyword whose value will be copied if the keyword is
!    not specified. Defaults to not set.
! At most one of default_value and default_keyword should be specified.
! is_optional specifies whether or not the keyword is required to be given.
!    is_optional defaults to .false. unless a default is given, when it
!    defaults to and is required to be .true..
! is_path specifies that the value will be a path to a file or directory.
!    Paths are converted to absolute paths (from /) automatically.
!    Defaults to .false..
! allowed_in_file specifies that the value may be read from or printed to file.
!    Defaults to .true..
! can_be_interactive specifies that the value may be set interactively by the
!    user.
!    Defaults to .true.
! flag specifies the flag by which the keyword can be alternately called.
!    Defaults to ' '.
function new_KeywordData(keyword,helptext,default_value,default_keyword, &
   & is_optional,is_path,allowed_in_file,can_be_interactive,             &
   & flag_without_arguments,flag_with_arguments) result(this)
  implicit none
  
  character(*), intent(in)           :: keyword
  character(*), intent(in)           :: helptext
  character(*), intent(in), optional :: default_value
  character(*), intent(in), optional :: default_keyword
  logical,      intent(in), optional :: is_optional
  logical,      intent(in), optional :: is_path
  logical,      intent(in), optional :: allowed_in_file
  logical,      intent(in), optional :: can_be_interactive
  character(1), intent(in), optional :: flag_without_arguments
  character(1), intent(in), optional :: flag_with_arguments
  type(KeywordData)                  :: this
  
  ! Check for incompatible arguments.
  if (present(is_optional)) then
    if (present(default_value) .or. present(default_keyword)) then
      call print_line(CODE_ERROR//': the argument "is_optional" should not be &
         &given for keywords which have defaults.')
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
  this%keyword = lower_case(keyword)
  this%helptext = helptext
  
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
  
  if (present(can_be_interactive)) then
    this%can_be_interactive = can_be_interactive
  else
    this%can_be_interactive = .true.
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
  else
    this%default_type_ = 0
  endif
  
  ! Unset value.
  call this%unset()
end function

! ----------------------------------------------------------------------
! Sets the keyword interactively, asking the user for a value.
! ----------------------------------------------------------------------
subroutine set_interactively(this)
  implicit none
  
  class(KeywordData), intent(inout) :: this
  
  type(String) :: input
  
  call print_line('')
  call print_line(this%helptext)
  if (this%is_set()) then
    call print_line(this%keyword//' currently has the value "'//this%value() &
       & //'".')
    call print_line('Please press <Enter> to accept this value, or enter a &
       & new value.')
    input = read_line_from_user()
    if (input/='') then
      call this%set(input)
    endif
  else
    if (this%default_type_==0) then
      call print_line(this%keyword//' is unset and has no default.')
      do while (.not. this%is_set())
        call print_line('Please enter a value.')
        input = read_line_from_user()
        if (input/='') then
          call this%set(input)
        endif
      enddo
    else
      if (this%default_type_==1) then
        call print_line(this%keyword//' is unset and is optional.')
        call print_line('Please press <Enter> to leave it unset, or enter a &
           &value.')
      elseif (this%default_type_==2) then
        call print_line(this%keyword//' is unset, and will default to the &
           &value "'//this%default_//'".')
        call print_line('Please press <Enter> to accept this value, or enter &
           &a value.')
      elseif (this%default_type_==3) then
        call print_line(this%keyword//' is unset, and will default to the &
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
end subroutine

! ----------------------------------------------------------------------
! Process and check value.
! ----------------------------------------------------------------------
! Throws and error if a value is required but has not been set.
! Turns paths into absolute form.
subroutine process_and_check(this)
  implicit none
  
  class(KeywordData), intent(inout) :: this
  
  if (this%default_type_==0 .and. .not. this%is_set()) then
    call print_line(ERROR//': the keyword '//this%keyword//' has not been &
       &set. this keyword is not optional.')
    stop
  endif
  
  if (this%is_path .and. this%is_set()) then
    this%value_ = format_path(this%value_)
  endif
end subroutine

! ----------------------------------------------------------------------
! Prints helptext and relevant keyword settings.
! ----------------------------------------------------------------------
subroutine print_help(this)
  implicit none
  
  class(KeywordData), intent(in) :: this
  
  call print_line('')
  call print_line(this%helptext)
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
end subroutine
end module

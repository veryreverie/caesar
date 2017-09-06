! ======================================================================
! A dictionary of keys and values, both of type String.
! ======================================================================
module dictionary_module
  use constants_module, only : dp
  use string_module
  use io_module
  use help_module
  implicit none
  
  private
  
  public :: size          ! The number of key/value pairs.
  public :: assignment(=) ! Make a Dictionary from an array of keywords.
  public :: operator(//)  ! Concatenate two Dictionaries.
  
  ! ------------------------------
  ! A dictionary of keys and values.
  ! ------------------------------
  type, public :: Dictionary
    type(KeywordData), allocatable :: keywords(:)
  contains
    ! ----------
    ! Getters.
    ! ----------
    ! Private function to return the index of a keyword.
    generic,   private :: index =>                    &
                        & index_Dictionary_character, &
                        & index_Dictionary_String
    procedure, private :: index_Dictionary_character
    procedure, private :: index_Dictionary_String
    
    ! As above, but by flag rather than by keyword.
    generic,   private :: index_by_flag =>                    &
                        & index_by_flag_Dictionary_character, &
                        & index_by_flag_Dictionary_String
    procedure, private :: index_by_flag_Dictionary_character
    procedure, private :: index_by_flag_Dictionary_String
    
    ! Returns the corresponding keyword.
    generic,   public  :: flag_to_keyword =>                    &
                        & flag_to_keyword_Dictionary_character, &
                        & flag_to_keyword_Dictionary_String
    procedure, private :: flag_to_keyword_Dictionary_character
    procedure, private :: flag_to_keyword_Dictionary_String
    
    ! Returns whether or not a keyword is set.
    generic,   public  :: is_set =>                    &
                        & is_set_Dictionary_character, &
                        & is_set_Dictionary_String
    procedure, private :: is_set_Dictionary_character
    procedure, private :: is_set_Dictionary_String
    
    ! Returns the value matching a keyword.
    generic,   public  :: value =>                    &
                        & value_Dictionary_character, &
                        & value_Dictionary_String
    procedure, private :: value_Dictionary_character
    procedure, private :: value_Dictionary_String
    
    ! ----------
    ! Setters.
    ! ----------
    ! Sets a keyword.
    generic,   public  :: set =>                    &
                        & set_Dictionary_character, &
                        & set_Dictionary_String
    procedure, private :: set_Dictionary_character
    procedure, private :: set_Dictionary_String
  
    ! Unsets a keyword.
    generic,   public  :: unset =>                    &
                        & unset_Dictionary_character, &
                        & unset_Dictionary_String
    procedure, private :: unset_Dictionary_character
    procedure, private :: unset_Dictionary_String
  
    ! Sets the value of a keyword.
    generic,   public  :: set_value =>                              &
                        & set_value_Dictionary_character_character, &
                        & set_value_Dictionary_character_String,    &
                        & set_value_Dictionary_String_character,    &
                        & set_value_Dictionary_String_String
    procedure, private :: set_value_Dictionary_character_character
    procedure, private :: set_value_Dictionary_character_String
    procedure, private :: set_value_Dictionary_String_character
    procedure, private :: set_value_Dictionary_String_String
  
    
    ! Appends to the value of a keyword.
    generic, public :: append_value =>                              &
                     & append_value_Dictionary_character_character, &
                     & append_value_Dictionary_character_String,    &
                     & append_value_Dictionary_String_character,    &
                     & append_value_Dictionary_String_String
    procedure, private :: append_value_Dictionary_character_character
    procedure, private :: append_value_Dictionary_character_String
    procedure, private :: append_value_Dictionary_String_character
    procedure, private :: append_value_Dictionary_String_String
    
    ! ----------
    ! File operations.
    ! ----------
    ! Write to a file.
    generic,   public  :: write_file =>                    &
                        & write_file_Dictionary_character, &
                        & write_file_Dictionary_String
    procedure, private :: write_file_Dictionary_character
    procedure, private :: write_file_Dictionary_String
  
    ! Read keywords from a file.
    generic,   public  :: read_file =>                    &
                        & read_file_Dictionary_character, &
                        & read_file_Dictionary_String
    procedure, private :: read_file_Dictionary_character
    procedure, private :: read_file_Dictionary_String
  
    
    ! ----------
    ! Set defaults, process paths and check non-optional keywords are set.
    ! ----------
    procedure, public  :: process_and_check_inputs => &
                        & process_and_check_inputs_Dictionary
    
  end type
  
  ! ------------------------------
  ! Procedures acting on a Dictionary.
  ! ------------------------------
  interface new
    module procedure new_Dictionary
  end interface
  
  interface size
    module procedure size_Dictionary
  end interface
  
  interface assignment(=)
    module procedure assign_Dictionary_KeywordDatas
  end interface
  
  interface operator(//)
    module procedure concatenate_Dictionary_Dictionary
  end interface
  
contains

! ----------------------------------------------------------------------
! Private allocate(Dictionary) subroutine.
! ----------------------------------------------------------------------
subroutine new_Dictionary(this,no_entries)
  implicit none
  
  type(Dictionary), intent(out) :: this
  integer,          intent(in)  :: no_entries
  
  integer :: ialloc
  
  allocate(this%keywords(no_entries), stat=ialloc); call err(ialloc)
end subroutine

! ----------------------------------------------------------------------
! size(Dictionary).
! ----------------------------------------------------------------------
function size_Dictionary(this) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: this
  integer                      :: output
  
  output = size(this%keywords)
end function

! ----------------------------------------------------------------------
! Takes an array of KeywordData and returns a Dictionary.
! ----------------------------------------------------------------------
subroutine assign_Dictionary_KeywordDatas(output,input)
  implicit none
  
  type(KeywordData), intent(in)  :: input(:)
  type(Dictionary),  intent(out) :: output
  
  output%keywords = input
end subroutine

! ----------------------------------------------------------------------
! Concatenates two Dictionaries.
! ----------------------------------------------------------------------
function concatenate_Dictionary_Dictionary(dict1,dict2) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: dict1
  type(Dictionary), intent(in) :: dict2
  type(Dictionary)             :: output
  
  call new(output,size(dict1)+size(dict2))
  output%keywords(:size(dict1)) = dict1%keywords
  output%keywords(size(dict1)+1:) = dict2%keywords
end function

! ----------------------------------------------------------------------
! Get the index where the keyword is stored.
! ----------------------------------------------------------------------
! Private function.
! Throws an error if the keyword is not found.
! If there are duplicate keys, returns the first match.
function index_Dictionary_character(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: keyword
  integer                       :: output
  
  integer :: i
  
  do i=1,size(this)
    if (this%keywords(i)%keyword == keyword) then
      output = i
      return
    endif
  enddo
  
  call print_line('Error: keyword '//keyword//' not found.')
  call err()
  
  ! Prevents a warning.
  output = 0
end function

function index_Dictionary_String(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: keyword
  integer                       :: output
  
  output = this%index(char(keyword))
end function

! ----------------------------------------------------------------------
! As above, but by flag rather than keyword.
! ----------------------------------------------------------------------
function index_by_flag_Dictionary_character(this,flag) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(1),      intent(in) :: flag
  integer                       :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this)
    if (this%keywords(i)%flag == flag) then
      output = i
      return
    endif
  enddo
  
  call print_line('Error: flag '//flag//' not found.')
  call err()
end function

function index_by_flag_Dictionary_String(this,flag) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: flag
  integer                       :: output
  
  output = this%index_by_flag(char(slice(flag,1,1)))
end function

! ----------------------------------------------------------------------
!  Returns the keyword corresponding to the given flag.
! ----------------------------------------------------------------------
function flag_to_keyword_Dictionary_character(this,flag) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(1),      intent(in) :: flag
  type(String)                  :: output
  
  output = this%keywords(this%index_by_flag(flag))%keyword
end function

function flag_to_keyword_Dictionary_String(this,flag) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: flag
  type(String)                  :: output
  
  output = this%flag_to_keyword(char(flag))
end function

! ----------------------------------------------------------------------
! Get whether or not a keyword has been set.
! ----------------------------------------------------------------------
function is_set_Dictionary_character(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: keyword
  logical                       :: output
  
  output = this%keywords(this%index(keyword))%is_set
end function

function is_set_Dictionary_String(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: keyword
  logical                       :: output
  
  output = this%is_set(char(keyword))
end function

! ----------------------------------------------------------------------
! Get the value corresponding to a given key.
! ----------------------------------------------------------------------
! Throws an error if the key has not been set or is boolean.
function value_Dictionary_character(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: keyword
  type(String)                  :: output
  
  integer :: i
  
  i = this%index(keyword)
  if (.not. this%keywords(i)%is_set) then
    call print_line('Error: keyword '//keyword//' has not been set.')
    call err()
  elseif (this%keywords(i)%is_boolean) then
    call print_line('Error: keyword '//keyword//' is boolean.')
    call err()
  endif
  
  output = this%keywords(i)%value
end function

function value_Dictionary_String(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: keyword
  type(String)                  :: output
  
  output = this%value(char(keyword))
end function

! ----------------------------------------------------------------------
! Sets a given keyword.
! ----------------------------------------------------------------------
subroutine set_Dictionary_character(this, keyword)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  
  this%keywords(this%index(keyword))%is_set = .true.
end subroutine

subroutine set_Dictionary_String(this, keyword)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(String),      intent(in)    :: keyword
  
  call this%set(char(keyword))
end subroutine

! ----------------------------------------------------------------------
! Unsets a given keyword.
! ----------------------------------------------------------------------
subroutine unset_Dictionary_character(this, keyword)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  
  integer :: i
  
  i = this%index(keyword)
  this%keywords(i)%keyword = ''
  this%keywords(i)%is_set = .false.
end subroutine

subroutine unset_Dictionary_String(this, keyword)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(String),      intent(in)    :: keyword
  
  call this%unset(char(keyword))
end subroutine

! ----------------------------------------------------------------------
! Sets the value corresponding to a given keyword.
! ----------------------------------------------------------------------
! If only_set_if_not_set is set, the value will not be overwritten if set.
subroutine set_value_Dictionary_character_character(this, keyword, value, &
   & only_set_if_not_set)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: keyword
  character(*),      intent(in)           :: value
  logical,           intent(in), optional :: only_set_if_not_set
  
  logical :: only_set ! = only_set_if_not_set.
  integer :: i
  
  if (present(only_set_if_not_set)) then
    only_set = only_set_if_not_set
  else
    only_set = .false.
  endif
  
  i = this%index(keyword)
  if (.not. (only_set .and. this%keywords(i)%is_set)) then
    this%keywords(i)%value = value
  endif
  this%keywords(i)%is_set = .true.
end subroutine

subroutine set_value_Dictionary_character_String(this, keyword, value, &
   & only_set_if_not_set)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: keyword
  type(String),      intent(in)           :: value
  logical,           intent(in), optional :: only_set_if_not_set
  
  if (present(only_set_if_not_set)) then
    call this%set_value(keyword, char(value), only_set_if_not_set)
  else
    call this%set_value(keyword, char(value))
  endif
end subroutine

subroutine set_value_Dictionary_String_character(this, keyword, value, &
   & only_set_if_not_set)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: keyword
  character(*),      intent(in)           :: value
  logical,           intent(in), optional :: only_set_if_not_set
  
  if (present(only_set_if_not_set)) then
    call this%set_value(char(keyword), value, only_set_if_not_set)
  else
    call this%set_value(char(keyword), value)
  endif
end subroutine

subroutine set_value_Dictionary_String_String(this, keyword, value, &
   & only_set_if_not_set)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: keyword
  type(String),      intent(in)           :: value
  logical,           intent(in), optional :: only_set_if_not_set
  
  if (present(only_set_if_not_set)) then
    call this%set_value(char(keyword), char(value), only_set_if_not_set)
  else
    call this%set_value(char(keyword), char(value))
  endif
end subroutine

! ----------------------------------------------------------------------
! Appends to the value corresponding to a given keyword.
! ----------------------------------------------------------------------
subroutine append_value_Dictionary_character_character(this, keyword, value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  character(*),      intent(in)    :: value
  
  integer :: i
  
  i = this%index(keyword)
  if (this%keywords(i)%is_set) then
    this%keywords(i)%value = this%keywords(i)%value // ' ' // value
  else
    this%keywords(i)%value = value
  endif
  this%keywords(i)%is_set = .true.
end subroutine

subroutine append_value_Dictionary_character_String(this, keyword, value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  type(String),      intent(in)    :: value
  
  call this%append_value(keyword, char(value))
end subroutine

subroutine append_value_Dictionary_String_character(this, keyword, value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(String),      intent(in)    :: keyword
  character(*),      intent(in)    :: value
  
  call this%append_value(char(keyword), value)
end subroutine

subroutine append_value_Dictionary_String_String(this, keyword, value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(String),      intent(in)    :: keyword
  type(String),      intent(in)    :: value
  
  call this%append_value(char(keyword), char(value))
end subroutine

! ----------------------------------------------------------------------
! Writes a Dictionary to file.
! ----------------------------------------------------------------------
subroutine write_file_Dictionary_character(this,filename)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: filename
  
  integer :: i
  integer :: dictionary_file
  
  dictionary_file = open_write_file(filename)
  do i=1,size(this)
    if (.not. this%keywords(i)%is_set) then
      cycle
    elseif (.not. this%keywords(i)%allowed_in_file) then
      cycle
    elseif (this%keywords(i)%is_boolean) then
      call print_line(dictionary_file, this%keywords(i)%keyword)
    else
      call print_line(dictionary_file, this%keywords(i)%keyword//' '// &
                                     & this%keywords(i)%value)
    endif
  enddo
  close(dictionary_file)
end subroutine

subroutine write_file_Dictionary_String(this,filename)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: filename
  
  call this%write_file(char(filename))
end subroutine

! ----------------------------------------------------------------------
! Reads Dictionary entries from a file.
! ----------------------------------------------------------------------
! The dictionary must already have been initialised from a list of keywords
!    before this is called.
! The contents of the file will be checked against the dictionary's keywords.
! If only_set_if_not_set is set, then only keywords which have not already
!    been set will be modified. This defaults to .false..
subroutine read_file_Dictionary_character(this, filename, &
   & only_set_if_not_set)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: filename
  logical,           intent(in), optional :: only_set_if_not_set
  
  ! Files.
  type(String), allocatable :: dictionary_file(:)
  
  ! Temporary variables.
  integer                   :: i,j
  type(String), allocatable :: line(:)
  logical                   :: only_set ! = only_set_if_not_set
  
  if (present(only_set_if_not_set)) then
    only_set = only_set_if_not_set
  else
    only_set = .false.
  endif
  
  ! Read file.
  dictionary_file = read_lines(filename)
  
  ! Process file.
  ! Each line is expected to be of the form '  key  value  ! comments  '
  do i=1,size(dictionary_file)
    ! Strip leading and trailing spaces.
    line = [trim(dictionary_file(i))]
    
    ! Ignore blank lines and comment lines (those beginning with a '!').
    if (len(line(1))==0) then
      cycle
    elseif (slice(line(1),1,1)=='!') then
      cycle
    endif
    
    ! Strip out anything after a comment character (!).
    line = split(line(1), '!')
    line = split(line(1))
    
    ! Find keyword in arguments.
    j = this%index(lower_case(line(1)))
    if (.not. this%keywords(j)%allowed_in_file) then
      call print_line('Error: the keyword '//this%keywords(j)%keyword// &
         & ' should not appear in input files.')
      call err()
    elseif (.not. (only_set .and. this%keywords(j)%is_set)) then
      if (this%keywords(j)%is_boolean) then
        if (size(line)>1) then
          call print_line('Error: the boolean keyword '// &
             & this%keywords(j)%keyword//' has been specified with an &
             &argument in file '//filename)
          call err()
        endif
      else
        if (size(line)==1) then
          this%keywords(j)%value = ''
        else
          this%keywords(j)%value = join(line(2:))
        endif
      endif
    endif
    this%keywords(j)%is_set = .true.
  enddo
end subroutine

subroutine read_file_Dictionary_String(this, filename, &
   & only_set_if_not_set)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: filename
  logical,           intent(in), optional :: only_set_if_not_set
  
  if (present(only_set_if_not_set)) then
    call this%read_file(char(filename), only_set_if_not_set)
  else
    call this%read_file(char(filename))
  endif
end subroutine

! ----------------------------------------------------------------------
! Set defaults for keywords which have not been set.
! Convert all paths to absolute format (from /).
! Check that all non-optional keywords have been set.
! ----------------------------------------------------------------------
subroutine process_and_check_inputs_Dictionary(this)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  
  integer :: i,j
  
  do_i : do i=1,size(this%keywords)
    if (.not. this%keywords(i)%is_set) then
      if (this%keywords(i)%default_value/='') then
        ! Set keywords with default values.
        this%keywords(i)%value = this%keywords(i)%default_value
        this%keywords(i)%is_set = .true.
      elseif (this%keywords(i)%default_keyword/='') then
        ! Set keywords which default to other keywords.
        do j=1,size(this%keywords)
          if (this%keywords(j)%keyword==this%keywords(i)%default_keyword) then
            if (.not. this%keywords(j)%is_set) then
              call print_line('Error: keyword '//this%keywords(i)%keyword // &
                 & 'defaults to keyword '//this%keywords(j)%keyword       // &
                 & '. This has not been set.')
              call err()
            endif
            this%keywords(i)%value = this%keywords(j)%value
            this%keywords(i)%is_set = .true.
            cycle do_i
          endif
        enddo
      elseif (.not. this%keywords(i)%is_optional) then
        ! Stop if a non-optional keyword has not been set.
        call print_line('Error: the non-optional keyword '// &
           & this%keywords(i)%keyword//' has not been set.')
        stop
      endif
    endif
    
    if (this%keywords(i)%is_path .and. this%keywords(i)%is_set) then
      ! Convert paths to absolute format (from /).
      this%keywords(i)%value = format_path(this%keywords(i)%value)
    endif
  enddo do_i
end subroutine
end module

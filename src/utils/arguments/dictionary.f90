! ======================================================================
! A dictionary of keys and values, both of type String.
! ======================================================================
module dictionary_submodule
  use precision_module
  use io_module
  use keyword_submodule
  use logic_module
  implicit none
  
  private
  
  public :: Dictionary
  public :: size ! The number of key/value pairs.
  
  ! ------------------------------
  ! A dictionary of keys and values.
  ! ------------------------------
  type :: Dictionary
    type(KeywordData), allocatable :: keywords(:)
  contains
    ! ----------
    ! Make a Dictionary from an array of keywords.
    ! ----------
    generic,   public  :: assignment(=) => assign_Dictionary_KeywordDatas
    procedure, private :: assign_Dictionary_KeywordDatas
    
    ! ----------
    ! Concatenate two Dictionaries.
    ! ----------
    generic,   public  :: operator(//) => concatenate_Dictionary_Dictionary
    procedure, private :: concatenate_Dictionary_Dictionary
    
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
    
    ! Returns whether or not a keyword is set with a value.
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
    ! Unsets a keyword.
    generic,   public  :: unset =>                    &
                        & unset_Dictionary_character, &
                        & unset_Dictionary_String
    procedure, private :: unset_Dictionary_character
    procedure, private :: unset_Dictionary_String
  
    ! Sets a keyword and sets its value.
    generic,   public  :: set =>                              &
                        & set_Dictionary_character_character, &
                        & set_Dictionary_character_String,    &
                        & set_Dictionary_String_character,    &
                        & set_Dictionary_String_String
    procedure, private :: set_Dictionary_character_character
    procedure, private :: set_Dictionary_character_String
    procedure, private :: set_Dictionary_String_character
    procedure, private :: set_Dictionary_String_String
  
    
    ! Appends to the value of a keyword.
    ! Returns an error if the keyword has no value set.
    generic, public :: append_to_value =>                              &
                     & append_to_value_Dictionary_character_character, &
                     & append_to_value_Dictionary_character_String,    &
                     & append_to_value_Dictionary_String_character,    &
                     & append_to_value_Dictionary_String_String
    procedure, private :: append_to_value_Dictionary_character_character
    procedure, private :: append_to_value_Dictionary_character_String
    procedure, private :: append_to_value_Dictionary_String_character
    procedure, private :: append_to_value_Dictionary_String_String
    
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
    procedure, public :: set_interactively => set_interactively_Dictionary
    procedure, public :: process_and_check_inputs => &
                       & process_and_check_inputs_Dictionary
    
  end type
  
  ! ------------------------------
  ! Procedures acting on a Dictionary.
  ! ------------------------------
  interface Dictionary
    module procedure new_Dictionary
  end interface
  
  interface size
    module procedure size_Dictionary
  end interface
contains

! ----------------------------------------------------------------------
! Private allocate(Dictionary) subroutine.
! ----------------------------------------------------------------------
function new_Dictionary(no_entries) result(this)
  implicit none
  
  integer, intent(in) :: no_entries
  type(Dictionary)    :: this
  
  integer :: ialloc
  
  allocate(this%keywords(no_entries), stat=ialloc); call err(ialloc)
end function

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
subroutine assign_Dictionary_KeywordDatas(this,that)
  implicit none
  
  class(Dictionary), intent(out) :: this
  type(KeywordData), intent(in)  :: that(:)
  
  this%keywords = that
end subroutine

! ----------------------------------------------------------------------
! Concatenates two Dictionaries.
! ----------------------------------------------------------------------
function concatenate_Dictionary_Dictionary(this,that) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  class(Dictionary), intent(in) :: that
  type(Dictionary)              :: output
  
  output%keywords = [this%keywords, that%keywords]
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
  
  output = first(this%keywords%keyword == keyword,default=0)
  
  if (output==0) then
    call print_line(ERROR//': unexpected keyword: '//keyword//'.')
    stop
  endif
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
  
  output = first(this%keywords,flag_matches,default=0)
  
  if (output==0) then
    call print_line(ERROR//': unexpected flag: '//flag//'.')
    stop
  endif
contains
  ! Lambda for checking if a flag matches.
  ! Captures flag from the function.
  function flag_matches(input) result(output)
    implicit none
    
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
! Get whether or not a keyword has been set with a value.
! ----------------------------------------------------------------------
function is_set_Dictionary_character(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: keyword
  logical                       :: output
  
  output = this%keywords(this%index(keyword))%is_set()
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
  
  output = this%keywords(this%index(keyword))%value()
end function

function value_Dictionary_String(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: keyword
  type(String)                  :: output
  
  output = this%value(char(keyword))
end function

! ----------------------------------------------------------------------
! Unsets a given keyword.
! ----------------------------------------------------------------------
subroutine unset_Dictionary_character(this, keyword)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  
  call this%keywords(this%index(keyword))%unset()
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
! If only_update_if_unset is set, the value will not be overwritten if set.
!    defaults to .false..
subroutine set_Dictionary_character_character(this,keyword,value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: keyword
  character(*),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  if (present(only_update_if_unset)) then
    call this%keywords(this%index(keyword))%set(value, only_update_if_unset) 
  else
    call this%keywords(this%index(keyword))%set(value)
  endif
end subroutine

subroutine set_Dictionary_character_String(this,keyword,value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: keyword
  type(String),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  if (present(only_update_if_unset)) then
    call this%set(keyword, char(value), only_update_if_unset)
  else
    call this%set(keyword, char(value))
  endif
end subroutine

subroutine set_Dictionary_String_character(this, keyword, value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: keyword
  character(*),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  if (present(only_update_if_unset)) then
    call this%set(char(keyword), value, only_update_if_unset)
  else
    call this%set(char(keyword), value)
  endif
end subroutine

subroutine set_Dictionary_String_String(this,keyword,value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: keyword
  type(String),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  if (present(only_update_if_unset)) then
    call this%set(char(keyword), char(value), only_update_if_unset)
  else
    call this%set(char(keyword), char(value))
  endif
end subroutine

! ----------------------------------------------------------------------
! Appends to the value corresponding to a given keyword.
! ----------------------------------------------------------------------
subroutine append_to_value_Dictionary_character_character(this,keyword,value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  character(*),      intent(in)    :: value
  
  call this%keywords(this%index(keyword))%append(value)
end subroutine

subroutine append_to_value_Dictionary_character_String(this,keyword,value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  type(String),      intent(in)    :: value
  
  call this%append_to_value(keyword, char(value))
end subroutine

subroutine append_to_value_Dictionary_String_character(this,keyword,value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(String),      intent(in)    :: keyword
  character(*),      intent(in)    :: value
  
  call this%append_to_value(char(keyword), value)
end subroutine

subroutine append_to_value_Dictionary_String_String(this,keyword,value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(String),      intent(in)    :: keyword
  type(String),      intent(in)    :: value
  
  call this%append_to_value(char(keyword), char(value))
end subroutine

! ----------------------------------------------------------------------
! Writes a Dictionary to file.
! ----------------------------------------------------------------------
subroutine write_file_Dictionary_character(this,filename)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: filename
  
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
    if (.not. this%keywords(i)%allowed_in_file) then
      cycle
    elseif (.not. this%keywords(i)%is_set()) then
      cycle
    else
      to_write = .true.
      max_length = max(max_length, len(this%keywords(i)%keyword))
    endif
  enddo
  
  if (to_write) then
    dictionary_file = filename
    do i=1,size(this)
      if (.not. this%keywords(i)%allowed_in_file) then
        cycle
      elseif (.not. this%keywords(i)%is_set()) then
        cycle
      else
        call dictionary_file%print_line(                           &
           & this%keywords(i)%keyword //                           &
           & spaces(max_length-len(this%keywords(i)%keyword)+1) // &
           & this%keywords(i)%value())
      endif
    enddo
  endif
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
! If only_update_if_unset is set, then only keywords which have not already
!    been set will be modified. This defaults to .false..
subroutine read_file_Dictionary_character(this, filename, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: filename
  logical,           intent(in), optional :: only_update_if_unset
  
  ! Files.
  type(IFile) :: dictionary_file
  
  ! Temporary variables.
  integer                   :: i,j
  type(String), allocatable :: line(:)
  logical                   :: only_if_unset ! = only_update_if_unset
  
  if (present(only_update_if_unset)) then
    only_if_unset = only_update_if_unset
  else
    only_if_unset = .false.
  endif
  
  ! Read file.
  dictionary_file = filename
  
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
    line = split(line(1), '!')
    line = split(line(1))
    
    ! Find keyword in arguments.
    j = this%index(lower_case(line(1)))
    if (.not. this%keywords(j)%allowed_in_file) then
      call print_line(ERROR//': the keyword '//this%keywords(j)%keyword// &
         & ' should not appear in input files.')
      call err()
    elseif (only_if_unset .and. this%keywords(j)%is_set()) then
      call print_line(WARNING//': the keyword '//this%keywords(j)%keyword// &
         & ' has been specified in multiple places.')
    else
      if (size(line)==1) then
        call print_line(ERROR//': the keyword '//this%keywords(j)%keyword// &
           & 'has been specified without a value.')
        stop
      else
        call this%keywords(j)%set(join(line(2:)),only_if_unset)
      endif
    endif
  enddo
end subroutine

subroutine read_file_Dictionary_String(this, filename, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: filename
  logical,           intent(in), optional :: only_update_if_unset
  
  if (present(only_update_if_unset)) then
    call this%read_file(char(filename), only_update_if_unset)
  else
    call this%read_file(char(filename))
  endif
end subroutine

! ----------------------------------------------------------------------
! Set defaults for keywords which have not been set.
! Convert all paths to absolute format (from /).
! Check that all non-optional keywords have been set.
! ----------------------------------------------------------------------
subroutine set_interactively_Dictionary(this)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  
  integer :: i
  
  do i=1,size(this)
    if (this%keywords(i)%can_be_interactive) then
      call this%keywords(i)%set_interactively()
    endif
  enddo
end subroutine

subroutine process_and_check_inputs_Dictionary(this)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  
  type(String) :: default_keyword
  
  integer :: i,j
  
  do_i : do i=1,size(this)
    if (.not. this%keywords(i)%is_set()) then
      default_keyword = this%keywords(i)%defaults_to_keyword()
      if (default_keyword == '') then
        call this%keywords(i)%set_default()
      else
        do j=1,size(this)
          if (this%keywords(j)%keyword==default_keyword) then
            if (this%keywords(j)%is_set()) then
              call this%keywords(i)%set(this%keywords(j)%value())
            endif
            cycle do_i
          endif
        enddo
        call err()
      endif
    endif
  enddo do_i
  
  do i=1,size(this)
    call this%keywords(i)%process_and_check()
  enddo
end subroutine
end module

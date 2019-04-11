! ======================================================================
! A dictionary of keys and values, both of type String.
! ======================================================================
module argument_dictionary_module
  use precision_module
  use abstract_module
  use io_module
  
  use keyword_module
  use common_keywords_module
  implicit none
  
  private
  
  public :: Dictionary
  public :: size
  public :: call_caesar
  
  ! ------------------------------
  ! A dictionary of keys and values.
  ! ------------------------------
  type, extends(NoDefaultConstructor) :: Dictionary
    ! Accepted keywords.
    type(KeywordData), allocatable, private :: keywords_(:)
  contains
    ! ----------
    ! Getters.
    ! ----------
    ! Private function to return the index of a keyword.
    generic,   private :: index =>                 &
                        & index_Dictionary_String, &
                        & index_Dictionary_character
    procedure, private :: index_Dictionary_String
    procedure, private :: index_Dictionary_character
    
    ! As above, but by flag rather than by keyword.
    generic,   private :: index_by_flag =>                    &
                        & index_by_flag_Dictionary_character, &
                        & index_by_flag_Dictionary_String
    procedure, private :: index_by_flag_Dictionary_character
    procedure, private :: index_by_flag_Dictionary_String
    
    ! Returns the corresponding keyword.
    generic,   private :: flag_to_keyword =>                    &
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
    
    ! Returns a list of python arguments.
    procedure, public :: python_arguments
    
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
                        & set_Dictionary_String_String,       &
                        & set_Dictionary_Dictionary
    procedure, private :: set_Dictionary_character_character
    procedure, private :: set_Dictionary_character_String
    procedure, private :: set_Dictionary_String_character
    procedure, private :: set_Dictionary_String_String
    procedure, private :: set_Dictionary_Dictionary
    
    ! Appends to the value of a keyword.
    ! Returns an error if the keyword has no value set.
    generic,   private :: append_to_value =>                              &
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
    procedure, public  :: set_interactively => set_interactively_Dictionary
    procedure, private :: process_and_check_inputs => &
                        & process_and_check_inputs_Dictionary
  end type
  
  interface Dictionary
    module procedure new_Dictionary_KeywordDatas
    module procedure new_Dictionary_arguments
  end interface
  
  interface size
    module procedure size_Dictionary
  end interface
  
  interface call_caesar
    module procedure call_caesar_String_Dictionary
    module procedure call_caesar_character_Dictionary
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
! Construct a Dictionary from and array of KeywordData.
function new_Dictionary_KeywordDatas(keywords) result(this)
  implicit none
  
  type(KeywordData), intent(in) :: keywords(:)
  type(Dictionary)              :: this
  
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
end function

! Construct a Dictionary from the command line input to Caesar.
function new_Dictionary_arguments(args,keywords_in) result(arguments)
  implicit none
  
  type(String),      intent(in) :: args(:)
  type(KeywordData), intent(in) :: keywords_in(:)
  type(Dictionary)              :: arguments
  
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
end function

! ----------------------------------------------------------------------
! size(Dictionary).
! ----------------------------------------------------------------------
function size_Dictionary(this) result(output)
  implicit none
  
  type(Dictionary), intent(in) :: this
  integer                      :: output
  
  output = size(this%keywords_)
end function

! ----------------------------------------------------------------------
! Get the index where the keyword is stored.
! ----------------------------------------------------------------------
! Private function.
! Throws an error if the keyword is not found.
! If there are duplicate keys, returns the first match.
impure elemental function index_Dictionary_String(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: keyword
  integer                       :: output
  
  output = first(this%keywords_%keyword() == keyword, default=0)
  
  if (output==0) then
    call print_line(ERROR//': unexpected keyword: '//keyword//'.')
    call quit()
  endif
end function

function index_Dictionary_character(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(*),      intent(in) :: keyword
  integer                       :: output
  
  output = this%index(str(keyword))
end function

! ----------------------------------------------------------------------
! As above, but by flag rather than keyword.
! ----------------------------------------------------------------------
function index_by_flag_Dictionary_character(this,flag) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  character(1),      intent(in) :: flag
  integer                       :: output
  
  output = first(this%keywords_,flag_matches,default=0)
  
  if (output==0) then
    call print_line(ERROR//': unexpected flag: '//flag//'.')
    call quit()
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

impure elemental function index_by_flag_Dictionary_String(this,flag) &
   & result(output)
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
  
  output = this%keywords_(this%index_by_flag(flag))%keyword()
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
  
  output = this%keywords_(this%index(keyword))%is_set()
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
  
  output = this%keywords_(this%index(keyword))%value()
end function

function value_Dictionary_String(this,keyword) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String),      intent(in) :: keyword
  type(String)                  :: output
  
  output = this%value(char(keyword))
end function

! ----------------------------------------------------------------------
! Return a list of python arguments.
! ----------------------------------------------------------------------
function python_arguments(this) result(output)
  implicit none
  
  class(Dictionary), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  integer :: i
  
  output = [String::]
  do i=1,size(this%keywords_)
    if ( this%keywords_(i)%is_set() .and.   &
       & this%keywords_(i)%pass_to_python() ) then
      output = [ output,                                                     &
               & this%keywords_(i)%keyword()//' '//this%keywords_(i)%value() ]
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Unsets a given keyword.
! ----------------------------------------------------------------------
subroutine unset_Dictionary_character(this, keyword)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  
  call this%keywords_(this%index(keyword))%unset()
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
  
  call this%keywords_(this%index(keyword))%set(value, only_update_if_unset) 
end subroutine

subroutine set_Dictionary_character_String(this,keyword,value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  character(*),      intent(in)           :: keyword
  type(String),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  call this%set(keyword, char(value), only_update_if_unset)
end subroutine

subroutine set_Dictionary_String_character(this, keyword, value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: keyword
  character(*),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  call this%set(char(keyword), value, only_update_if_unset)
end subroutine

subroutine set_Dictionary_String_String(this,keyword,value, &
   & only_update_if_unset)
  implicit none
  
  class(Dictionary), intent(inout)        :: this
  type(String),      intent(in)           :: keyword
  type(String),      intent(in)           :: value
  logical,           intent(in), optional :: only_update_if_unset
  
  call this%set(char(keyword), char(value), only_update_if_unset)
end subroutine

! ----------------------------------------------------------------------
! Sets all values from another dictionary.
! ----------------------------------------------------------------------
subroutine set_Dictionary_Dictionary(this,that)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  type(Dictionary),  intent(in)    :: that
  
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
end subroutine

! ----------------------------------------------------------------------
! Appends to the value corresponding to a given keyword.
! ----------------------------------------------------------------------
subroutine append_to_value_Dictionary_character_character(this,keyword,value)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  character(*),      intent(in)    :: keyword
  character(*),      intent(in)    :: value
  
  call this%keywords_(this%index(keyword))%append(value)
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
        call this%keywords_(j)%set(join(line(2:)),only_if_unset)
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
  
  call this%read_file(char(filename), only_update_if_unset)
end subroutine

! ----------------------------------------------------------------------
! Set defaults for keywords which have not been set.
! Convert all paths to absolute format (from /).
! Check that all non-optional keywords have been set.
! ----------------------------------------------------------------------
subroutine set_interactively_Dictionary(this)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  
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
end subroutine

subroutine process_and_check_inputs_Dictionary(this)
  implicit none
  
  class(Dictionary), intent(inout) :: this
  
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
end subroutine

! ----------------------------------------------------------------------
! Calls Caesar with a given dictionary of arguments.
! ----------------------------------------------------------------------
subroutine call_caesar_String_Dictionary(mode,arguments)
  implicit none
  
  type(String),     intent(in) :: mode
  type(Dictionary), intent(in) :: arguments
  
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
end subroutine

subroutine call_caesar_character_Dictionary(mode,arguments)
  implicit none
  
  character(*),     intent(in) :: mode
  type(Dictionary), intent(in) :: arguments
  
  call call_caesar(str(mode), arguments)
end subroutine
end module

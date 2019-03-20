! ======================================================================
! A dictionary of modes, containing:
!    - A brief description of each mode.
!    - The keywords for each mode.
!    - The helptext for said keyowrds.
! ======================================================================
module caesar_modes_module
  use io_module
  use abstract_module
  
  use keyword_module
  use argument_dictionary_module
  implicit none
  
  private
  
  public :: CaesarMode
  public :: CaesarModes
  public :: MainSubroutine
  public :: Dictionary
  public :: assignment(=)
  public :: add_mode
  
  ! An interface for the main subroutines of Caesar, each of which takes a
  !    dictionary of arguments and returns nothing.
  abstract interface
    subroutine MainSubroutine(arguments)
      import Dictionary
      implicit none
      
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  ! A container for a mode.
  type :: CaesarMode
    type(String)                               :: mode_name
    type(String)                               :: description
    type(KeywordData),         allocatable     :: keywords(:)
    procedure(MainSubroutine), pointer, nopass :: main_subroutine => null ()
    
    logical :: suppress_from_helptext = .false.
    logical :: suppress_settings_file = .false.
  contains
    procedure :: print_help => print_help_CaesarMode
  end type
  
  interface assignment(=)
    module procedure :: assign_CaesarMode
  end interface
  
  interface CaesarMode
    module procedure new_CaesarMode_character
    module procedure new_CaesarMode_String
    module procedure new_CaesarMode_character_character
    module procedure new_CaesarMode_character_String
    module procedure new_CaesarMode_String_character
    module procedure new_CaesarMode_String_String
  end interface
  
  ! A dictionary of modes.
  type :: CaesarModes
    type(CaesarMode), private, allocatable :: modes_(:)
  contains
    generic,   public  :: mode =>      &
                        & mode_String, &
                        & mode_character
    procedure, private :: mode_String
    procedure, private :: mode_character
    
    procedure, public  :: print_help => print_help_CaesarModes
  end type
  
  interface CaesarModes
    module procedure new_CaesarModes
    module procedure new_CaesarModes_CaesarMode
  end interface
  
  interface Dictionary
    module procedure new_Dictionary_CaesarMode
  end interface
  
  ! An array of modes, which will be populated at startup.
  type(CaesarModes) :: CAESAR_MODES
  
  interface add_mode
    module procedure add_mode_CaesarMode
  end interface
contains

! ----------------------------------------------------------------------
! CaesarMode procedures.
! ----------------------------------------------------------------------

! Constructor which picks the requested mode from CAESAR_MODES.
function new_CaesarMode_character(mode) result(this)
  implicit none
  
  character(*), intent(in) :: mode
  type(CaesarMode)         :: this
  
  this = CAESAR_MODES%mode(mode)
end function

function new_CaesarMode_String(mode) result(this)
  implicit none
  
  type(String), intent(in) :: mode
  type(CaesarMode)         :: this
  
  this = CaesarMode(char(mode))
end function

! Constructor for CaesarMode type. Takes:
!    - The name of the mode, e.g. 'harmonic'. This is converted to lower case.
!    - A brief description of the mode.
!    - The keywords associated with the mode.
!    - A pointer to the mode's subroutine.
! Can have suppress_from_helptext set, which stops the mode from appearing
!    in non-mode-specific help. This defaults to false.
function new_CaesarMode_character_character(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext,suppress_settings_file)      &
   & result(output)
  implicit none
  
  character(*),              intent(in)           :: mode_name
  character(*),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  logical,                   intent(in), optional :: suppress_settings_file
  type(CaesarMode)                                :: output
  
  output%mode_name = lower_case(mode_name)
  output%description = description
  output%keywords = keywords
  output%main_subroutine => main_subroutine
  
  if (present(suppress_from_helptext)) then
    output%suppress_from_helptext = suppress_from_helptext
  else
    output%suppress_from_helptext = .false.
  endif
  
  if (present(suppress_settings_file)) then
    output%suppress_settings_file = suppress_settings_file
  else
    output%suppress_settings_file = .false.
  endif
end function

function new_CaesarMode_character_String(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext,suppress_settings_file)      &
   & result(output)
  implicit none
  
  character(*),              intent(in)           :: mode_name
  type(String),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  logical,                   intent(in), optional :: suppress_settings_file
  type(CaesarMode)                                :: output
  
  output = CaesarMode( mode_name,              &
                     & char(description),      &
                     & keywords,               &
                     & main_subroutine,        &
                     & suppress_from_helptext, &
                     & suppress_settings_file)
end function

function new_CaesarMode_String_character(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext,suppress_settings_file)      &
   & result(output)
  implicit none
  
  type(String),              intent(in)           :: mode_name
  character(*),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  logical,                   intent(in), optional :: suppress_settings_file
  type(CaesarMode)                                :: output
  
  output = CaesarMode( char(mode_name),        &
                     & description,            &
                     & keywords,               &
                     & main_subroutine,        &
                     & suppress_from_helptext, &
                     & suppress_settings_file)
end function

function new_CaesarMode_String_String(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext,suppress_settings_file)      &
   & result(output)
  implicit none
  
  type(String),              intent(in)           :: mode_name
  type(String),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  logical,                   intent(in), optional :: suppress_settings_file
  type(CaesarMode)                                :: output
  
  output = CaesarMode( char(mode_name),        &
                     & char(description),      &
                     & keywords,               &
                     & main_subroutine,        &
                     & suppress_from_helptext, &
                     & suppress_settings_file)
end function

! Construct a dictionary from a CaesarMode.
function new_Dictionary_CaesarMode(mode) result(this)
  implicit none
  
  type(CaesarMode), intent(in) :: mode
  type(Dictionary)             :: this
  
  this = Dictionary(mode%keywords)
end function

! Copies one CaesarMode to another.
! Required to be explicit due to an apparent nagfort 6.2 bug with automatic
!    reallocation from an empty list.
! N.B. output%keywords = [KeywordData::] fails under nagfort 6.2.
subroutine assign_CaesarMode(output,input)
  implicit none
  
  class(CaesarMode), intent(out) :: output
  class(CaesarMode), intent(in)  :: input
  
  integer :: ialloc
  
  output%mode_name = input%mode_name
  output%description = input%description
  if (size(input%keywords)==0) then
    allocate(output%keywords(0), stat=ialloc); call err(ialloc)
  else
    output%keywords = input%keywords
  endif
  output%main_subroutine => input%main_subroutine
  output%suppress_from_helptext = input%suppress_from_helptext
  output%suppress_settings_file = input%suppress_settings_file
end subroutine

! Prints helptext.
subroutine print_help_CaesarMode(this)
  implicit none
  
  class(CaesarMode), intent(in) :: this
  
  call print_line('')
  call print_line(colour(this%mode_name,'cyan'))
  call print_line(this%description)
end subroutine

! ----------------------------------------------------------------------
! CaesarModes procedures.
! ----------------------------------------------------------------------

! Constructor for CaesarModes.
function new_CaesarModes() result(this)
  implicit none
  
  type(CaesarModes) :: this
  
  if (.not. allocated(CAESAR_MODES%modes_)) then
    call print_line(CODE_ERROR//': CAESAR_MODES has not been allocated.')
    call err()
  else
    this = CAESAR_MODES
  endif
end function

function new_CaesarModes_CaesarMode(modes) result(this)
  implicit none
  
  type(CaesarMode), intent(in) :: modes(:)
  type(CaesarModes)            :: this
  
  this%modes_ = modes
end function

! Returns the mode with a given name.
function mode_String(this,mode_name) result(output)
  implicit none
  
  class(CaesarModes), intent(in) :: this
  type(String),       intent(in) :: mode_name
  type(CaesarMode)               :: output
  
  type(String) :: lower_case_mode_name
  
  integer :: i
  
  lower_case_mode_name = lower_case(mode_name)
  
  i = first(this%modes_%mode_name==lower_case_mode_name, default=0)
  
  if (i==0) then
    call print_line(ERROR//': Unrecognised mode: '//mode_name)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  endif
  
  output = this%modes_(i)
end function

function mode_character(this,mode_name) result(output)
  implicit none
  
  class(CaesarModes), intent(in) :: this
  character(*),       intent(in) :: mode_name
  type(CaesarMode)               :: output
  
  output = this%mode(str(mode_name))
end function

! Prints helptext.
subroutine print_help_CaesarModes(this)
  implicit none
  
  class(CaesarModes), intent(in) :: this
  
  integer :: i
  
  do i=1,size(this%modes_)
    call this%modes_(i)%print_help()
  enddo
end subroutine

! Add a mode to the list of modes.
subroutine add_mode_CaesarMode(mode)
  implicit none
  
  type(CaesarMode), intent(in) :: mode
  
  if (.not. allocated(CAESAR_MODES%modes_)) then
    CAESAR_MODES%modes_ = [mode]
  else
    if (any(CAESAR_MODES%modes_%mode_name==mode%mode_name)) then
      call print_line(CODE_ERROR//': Trying to add the same mode twice.')
      call err()
    else
      CAESAR_MODES%modes_ = [CAESAR_MODES%modes_, mode]
    endif
  endif
end subroutine
end module

! ======================================================================
! A dictionary of modes, containing:
!    - A brief description of each mode.
!    - The keywords for each mode.
!    - The helptext for said keyowrds.
! ======================================================================
module caesar_modes_module
  use string_module
  use io_module
  
  use keyword_module
  implicit none
  
  ! An interface for the main subroutines of Caesar, each of which takes a
  !    dictionary of arguments and returns nothing.
  abstract interface
    subroutine MainSubroutine(arguments)
      use dictionary_module
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
  contains
    procedure :: print_help => print_help_CaesarMode
  end type
  
  interface CaesarMode
    module procedure new_CaesarMode_character_character
    module procedure new_CaesarMode_character_String
    module procedure new_CaesarMode_String_character
    module procedure new_CaesarMode_String_String
  end interface
  
  ! A dictionary of modes.
  type :: CaesarModes
    type(CaesarMode), private, allocatable :: modes_(:)
  contains
    generic,   public  :: mode => mode_character, &
                                & mode_String
    procedure, private ::         mode_character
    procedure, private ::         mode_String
    
    procedure, public  :: print_help => print_help_CaesarModes
  end type
  
  interface CaesarModes
    module procedure new_CaesarModes
  end interface
contains

! ----------------------------------------------------------------------
! CaesarMode procedures.
! ----------------------------------------------------------------------

! Constructor for CaesarMode type. Takes:
!    - The name of the mode, e.g. 'harmonic'. This is converted to lower case.
!    - A brief description of the mode.
!    - The keywords associated with the mode.
!    - A pointer to the mode's subroutine.
! Can have suppress_from_helptext set, which stops the mode from appearing
!    in non-mode-specific help. This defaults to false.
function new_CaesarMode_character_character(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext) result(output)
  implicit none
  
  character(*),              intent(in)           :: mode_name
  character(*),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
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
end function

function new_CaesarMode_character_String(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext) result(output)
  implicit none
  
  character(*),              intent(in)           :: mode_name
  type(String),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  type(CaesarMode)                                :: output
  
  if (present(suppress_from_helptext)) then
    output = CaesarMode( mode_name,         &
                       & char(description), &
                       & keywords,          &
                       & main_subroutine,   &
                       & suppress_from_helptext)
  else
    output = CaesarMode( mode_name,         &
                       & char(description), &
                       & keywords,          &
                       & main_subroutine)
  endif
end function

function new_CaesarMode_String_character(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext) result(output)
  implicit none
  
  type(String),              intent(in)           :: mode_name
  character(*),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  type(CaesarMode)                                :: output
  
  if (present(suppress_from_helptext)) then
    output = CaesarMode( char(mode_name), &
                       & description,     &
                       & keywords,        &
                       & main_subroutine, &
                       & suppress_from_helptext)
  else
    output = CaesarMode( char(mode_name), &
                       & description,     &
                       & keywords,        &
                       & main_subroutine)
  endif
end function

function new_CaesarMode_String_String(mode_name,description,keywords, &
   & main_subroutine,suppress_from_helptext) result(output)
  implicit none
  
  type(String),              intent(in)           :: mode_name
  type(String),              intent(in)           :: description
  type(KeywordData),         intent(in)           :: keywords(:)
  procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
  logical,                   intent(in), optional :: suppress_from_helptext
  type(CaesarMode)                                :: output
  
  if (present(suppress_from_helptext)) then
    output = CaesarMode( char(mode_name),   &
                       & char(description), &
                       & keywords,          &
                       & main_subroutine,   &
                       & suppress_from_helptext)
  else
    output = CaesarMode( char(mode_name),   &
                       & char(description), &
                       & keywords,          &
                       & main_subroutine)
  endif
end function

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
function new_CaesarModes(modes) result(output)
  implicit none
  
  type(CaesarMode), intent(in) :: modes(:)
  type(CaesarModes)            :: output
  
  output%modes_ = modes
end function

! Returns the mode with a given name.
function mode_character(this,mode_name) result(output)
  implicit none
  
  class(CaesarModes), intent(in) :: this
  character(*),       intent(in) :: mode_name
  type(CaesarMode)               :: output
  
  type(String) :: lower_case_name
  logical      :: success
  
  integer :: i
  
  success = .false.
  lower_case_name = lower_case(mode_name)
  do i=1,size(this%modes_)
    if (this%modes_(i)%mode_name==lower_case_name) then
      output = this%modes_(i)
      success = .true.
      exit
    endif
  enddo
  
  if (.not. success) then
    call print_line(colour('Error: unrecognised mode: ','red')//mode_name)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop
  endif
end function

function mode_String(this,mode_name) result(output)
  implicit none
  
  class(CaesarModes), intent(in) :: this
  type(String),       intent(in) :: mode_name
  type(CaesarMode)               :: output
  
  output = this%mode(char(mode_name))
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
end module

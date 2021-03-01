! ======================================================================
! A dictionary of modes, containing:
!    - A brief description of each mode.
!    - The keywords for each mode.
!    - The helptext for said keyowrds.
! ======================================================================
module caesar_caesar_modes_module
  use caesar_io_module
  use caesar_abstract_module
  
  use caesar_keyword_module
  use caesar_dictionary_module
  implicit none
  
  private
  
  public :: CaesarMode
  public :: CaesarModes
  public :: MainSubroutine
  public :: Dictionary
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
    procedure, public :: print_help => print_help_CaesarMode
    
    generic,   public  :: remove_keyword =>      &
                        & remove_keyword_String, &
                        & remove_keyword_character
    procedure, private :: remove_keyword_String
    procedure, private :: remove_keyword_character
  end type
  
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
  
  interface CaesarMode
    ! ----------------------------------------------------------------------
    ! CaesarMode procedures.
    ! ----------------------------------------------------------------------
    
    ! Constructor which picks the requested mode from CAESAR_MODES.
    module function new_CaesarMode_character(mode) result(this) 
      character(*), intent(in) :: mode
      type(CaesarMode)         :: this
    end function
  
    module function new_CaesarMode_String(mode) result(this) 
      type(String), intent(in) :: mode
      type(CaesarMode)         :: this
    end function
  
    ! Constructor for CaesarMode type. Takes:
    !    - The name of the mode, e.g. 'harmonic'. This is converted to lower case.
    !    - A brief description of the mode.
    !    - The keywords associated with the mode.
    !    - A pointer to the mode's module subroutine.
    ! Can have suppress_from_helptext set, which stops the mode from appearing
    !    in non-mode-specific help. This defaults to false.
    module function new_CaesarMode_character_character(mode_name,     &
        description,keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      character(*),              intent(in)           :: mode_name
      character(*),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(CaesarMode)                                :: output
    end function
  
    module function new_CaesarMode_character_String(mode_name,description, &
        keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      character(*),              intent(in)           :: mode_name
      type(String),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(CaesarMode)                                :: output
    end function
  
    module function new_CaesarMode_String_character(mode_name,description, &
        keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      type(String),              intent(in)           :: mode_name
      character(*),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(CaesarMode)                                :: output
    end function
  
    module function new_CaesarMode_String_String(mode_name,description, &
        keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      type(String),              intent(in)           :: mode_name
      type(String),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(CaesarMode)                                :: output
    end function
  end interface
  
  interface
    ! Remove a keyword from a CaesarMode.
    impure elemental module subroutine remove_keyword_String(this,keyword) 
      class(CaesarMode), intent(inout) :: this
      type(string),      intent(in)    :: keyword
    end subroutine
  end interface
  
  interface
    module subroutine remove_keyword_character(this,keyword) 
      class(CaesarMode), intent(inout) :: this
      character(*),      intent(in)    :: keyword
    end subroutine
  end interface
  
  interface Dictionary
    ! Construct a dictionary from a CaesarMode.
    module function new_Dictionary_CaesarMode(mode) result(this) 
      type(CaesarMode), intent(in) :: mode
      type(Dictionary)             :: this
    end function
  end interface
  
  interface print_help
    ! Prints helptext.
    module subroutine print_help_CaesarMode(this) 
      class(CaesarMode), intent(in) :: this
    end subroutine
  end interface
  
  interface CaesarModes
    ! ----------------------------------------------------------------------
    ! CaesarModes procedures.
    ! ----------------------------------------------------------------------
    
    ! Constructor for CaesarModes.
    module function new_CaesarModes() result(this) 
      type(CaesarModes) :: this
    end function
  
    module function new_CaesarModes_CaesarMode(modes) result(this) 
      type(CaesarMode), intent(in) :: modes(:)
      type(CaesarModes)            :: this
    end function
  end interface
  
  interface
    ! Returns the mode with a given name.
    module function mode_String(this,mode_name) result(output) 
      class(CaesarModes), intent(in) :: this
      type(String),       intent(in) :: mode_name
      type(CaesarMode)               :: output
    end function
  end interface
  
  interface
    module function mode_character(this,mode_name) result(output) 
      class(CaesarModes), intent(in) :: this
      character(*),       intent(in) :: mode_name
      type(CaesarMode)               :: output
    end function
  end interface
  
  interface
    ! Prints helptext.
    module subroutine print_help_CaesarModes(this) 
      class(CaesarModes), intent(in) :: this
    end subroutine
  end interface
  
  interface add_mode
    ! Add a mode to the list of modes.
    module subroutine add_mode_CaesarMode(mode) 
      type(CaesarMode), intent(in) :: mode
    end subroutine
  end interface
end module

! ======================================================================
! A dictionary of modes, containing:
!    - A brief description of each mode.
!    - The keywords for each mode.
!    - The helptext for said keyowrds.
! ======================================================================
module caesar_program_mode_module
  use caesar_io_module
  use caesar_abstract_module
  
  use caesar_keyword_module
  use caesar_dictionary_module
  implicit none
  
  private
  
  public :: ProgramMode
  public :: MainSubroutine
  public :: Dictionary
  
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
  type :: ProgramMode
    type(String)                               :: mode_name
    type(String)                               :: description
    type(KeywordData),         allocatable     :: keywords(:)
    procedure(MainSubroutine), pointer, nopass :: main_subroutine => null ()
    
    logical :: suppress_from_helptext = .false.
    logical :: suppress_settings_file = .false.
  contains
    procedure, public :: print_help => print_help_ProgramMode
    
    generic,   public  :: remove_keyword =>      &
                        & remove_keyword_String, &
                        & remove_keyword_character
    procedure, private :: remove_keyword_String
    procedure, private :: remove_keyword_character
  end type
  
  interface ProgramMode
    ! ----------------------------------------------------------------------
    ! ProgramMode procedures.
    ! ----------------------------------------------------------------------
  
    ! Constructor for ProgramMode type. Takes:
    !    - The name of the mode, e.g. 'harmonic'. This is converted to lower case.
    !    - A brief description of the mode.
    !    - The keywords associated with the mode.
    !    - A pointer to the mode's module subroutine.
    ! Can have suppress_from_helptext set, which stops the mode from appearing
    !    in non-mode-specific help. This defaults to false.
    module function new_ProgramMode_character_character(mode_name,   &
        description,keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      character(*),              intent(in)           :: mode_name
      character(*),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(ProgramMode)                               :: output
    end function
  
    module function new_ProgramMode_character_String(mode_name,description, &
        keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      character(*),              intent(in)           :: mode_name
      type(String),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(ProgramMode)                               :: output
    end function
  
    module function new_ProgramMode_String_character(mode_name,description, &
        keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      type(String),              intent(in)           :: mode_name
      character(*),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(ProgramMode)                               :: output
    end function
  
    module function new_ProgramMode_String_String(mode_name,description, &
        keywords,main_subroutine,suppress_from_helptext, &
          & suppress_settings_file) result(output) 
      type(String),              intent(in)           :: mode_name
      type(String),              intent(in)           :: description
      type(KeywordData),         intent(in)           :: keywords(:)
      procedure(MainSubroutine), intent(in), pointer  :: main_subroutine
      logical,                   intent(in), optional :: suppress_from_helptext
      logical,                   intent(in), optional :: suppress_settings_file
      type(ProgramMode)                               :: output
    end function
  end interface
  
  interface
    ! Remove a keyword from a ProgramMode.
    impure elemental module subroutine remove_keyword_String(this,keyword) 
      class(ProgramMode), intent(inout) :: this
      type(string),       intent(in)    :: keyword
    end subroutine
  end interface
  
  interface
    module subroutine remove_keyword_character(this,keyword) 
      class(ProgramMode), intent(inout) :: this
      character(*),       intent(in)    :: keyword
    end subroutine
  end interface
  
  interface Dictionary
    ! Construct a dictionary from a ProgramMode.
    module function new_Dictionary_ProgramMode(mode) result(this) 
      type(ProgramMode), intent(in) :: mode
      type(Dictionary)              :: this
    end function
  end interface
  
  interface print_help
    ! Prints helptext.
    module subroutine print_help_ProgramMode(this) 
      class(ProgramMode), intent(in) :: this
    end subroutine
  end interface
end module

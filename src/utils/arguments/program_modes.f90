module caesar_program_modes_module
  use caesar_io_module
  use caesar_abstract_module
  
  use caesar_keyword_module
  use caesar_dictionary_module
  use caesar_program_mode_module
  implicit none
  
  private
  
  public :: ProgramModes
  public :: add_mode
  
  ! A dictionary of modes.
  type :: ProgramModes
    type(ProgramMode), private, allocatable :: modes_(:)
  contains
    generic,   public  :: mode =>      &
                        & mode_String, &
                        & mode_character
    procedure, private :: mode_String
    procedure, private :: mode_character
    
    procedure, public  :: print_help => print_help_ProgramModes
  end type
  
  interface ProgramModes
    ! ----------------------------------------------------------------------
    ! ProgramModes procedures.
    ! ----------------------------------------------------------------------
    
    ! Constructor for ProgramModes.
    module function new_ProgramModes_ProgramMode(modes) result(this) 
      type(ProgramMode), intent(in) :: modes(:)
      type(ProgramModes)            :: this
    end function
  end interface
  
  interface
    ! Prints helptext.
    module subroutine print_help_ProgramModes(this) 
      class(ProgramModes), intent(in) :: this
    end subroutine
  end interface
  
  interface add_mode
    ! Add a mode to the list of modes.
    module subroutine add_mode_ProgramMode(mode) 
      type(ProgramMode), intent(in) :: mode
    end subroutine
  end interface
  
  interface
    ! Returns the mode with a given name.
    module function mode_String(this,mode_name) result(output) 
      class(ProgramModes), intent(in) :: this
      type(String),        intent(in) :: mode_name
      type(ProgramMode)               :: output
    end function
  end interface
  
  interface
    module function mode_character(this,mode_name) result(output) 
      class(ProgramModes), intent(in) :: this
      character(*),        intent(in) :: mode_name
      type(ProgramMode)               :: output
    end function
  end interface
end module

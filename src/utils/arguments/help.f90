! ======================================================================
! Handles input keywords and helptext.
! ======================================================================
module caesar_help_module
  use caesar_foundations_module
  use caesar_io_module
  use caesar_keyword_module
  use caesar_program_mode_module
  use caesar_program_modes_module
  implicit none
  
  private
  
  ! Prints help text.
  public :: help
  
  interface help
    ! ----------------------------------------------------------------------
    ! Prints helptext.
    ! ----------------------------------------------------------------------
      
    ! Print default helptext.
    module subroutine help_default(program_modes) 
      type(ProgramModes), intent(in) :: program_modes
    end subroutine
  
    ! Prints the helptext for a particular mode or keyword.
    module subroutine help_specific(keyword,program_mode) 
      type(String),      intent(in) :: keyword
      type(ProgramMode), intent(in) :: program_mode
    end subroutine
  end interface
end module

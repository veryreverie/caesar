! ======================================================================
! Handles input keywords and helptext.
! ======================================================================
module caesar_help_module
  use caesar_foundations_module
  use caesar_io_module
  use caesar_keyword_module
  use caesar_caesar_modes_module
  implicit none
  
  private
  
  ! Prints help text.
  public :: help
  
  interface help
    ! ----------------------------------------------------------------------
    ! Prints helptext.
    ! ----------------------------------------------------------------------
      
    ! Print default helptext.
    module subroutine help_default() 
    end subroutine
  
    ! Prints the helptext for a particular mode or keyword.
    module subroutine help_specific(keyword,mode) 
      type(String), intent(in) :: keyword
      type(String), intent(in) :: mode
    end subroutine
  end interface
end module

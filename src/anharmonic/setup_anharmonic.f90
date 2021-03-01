! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module caesar_setup_anharmonic_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  implicit none
  
  private
  
  public :: startup_setup_anharmonic
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_setup_anharmonic() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The main program.
    ! ----------------------------------------------------------------------
    module subroutine setup_anharmonic_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

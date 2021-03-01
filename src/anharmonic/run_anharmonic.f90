! ======================================================================
! Runs anharmonic DFT calculations, as set up by setup_anharmonic.
! ======================================================================
module caesar_run_anharmonic_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  implicit none
  
  private
  
  public :: startup_run_anharmonic
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_run_anharmonic() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine run_anharmonic_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

! ======================================================================
! The second stage of Caesar.
! Runs DFT for harmonic calculations.
! ======================================================================
module caesar_run_harmonic_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_run_harmonic
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_run_harmonic() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine run_harmonic_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

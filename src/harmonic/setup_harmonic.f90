! ======================================================================
! The first stage of Caesar.
! Generates supercells, and prepares harmonic DFT calculations.
! ======================================================================
module caesar_setup_harmonic_module
  use caesar_common_module
  
  use caesar_harmonic_data_module
  use caesar_generate_supercells_module
  implicit none
  
  private
  
  public :: startup_setup_harmonic
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_setup_harmonic() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine setup_harmonic_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

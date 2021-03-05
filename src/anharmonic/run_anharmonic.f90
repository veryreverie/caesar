! ======================================================================
! Runs anharmonic DFT calculations, as set up by setup_anharmonic.
! ======================================================================
module caesar_run_anharmonic_module
  use caesar_common_module
  use caesar_harmonic_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  implicit none
  
  private
  
  public :: run_anharmonic_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function run_anharmonic_mode() result(output)
      type(ProgramMode) :: output
    end function
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

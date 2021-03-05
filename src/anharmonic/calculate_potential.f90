! ======================================================================
! Calculates the anharmonic potential, using the results of run_anharmonic.
! ======================================================================
module caesar_calculate_potential_module
  use caesar_common_module
  use caesar_harmonic_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  use caesar_interpolation_module
  implicit none
  
  private
  
  public :: calculate_potential_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function calculate_potential_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine calculate_potential_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

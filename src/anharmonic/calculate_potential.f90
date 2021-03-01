! ======================================================================
! Calculates the anharmonic potential, using the results of run_anharmonic.
! ======================================================================
module caesar_calculate_potential_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  use caesar_interpolation_module
  implicit none
  
  private
  
  public :: startup_calculate_potential
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_calculate_potential() 
    end subroutine
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

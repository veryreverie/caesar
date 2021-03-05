! ======================================================================
! Calculates, under the harmonic approximation:
!    - The phonon dispersion curve along a specified path.
!    - The phonon density of states.
!    - The energy, free energy and entropy per unit cell.
! ======================================================================
! Should be run after calculate_normal_modes.
! Not required for anharmonic calculations.
module caesar_calculate_harmonic_observables_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: calculate_harmonic_observables_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function calculate_harmonic_observables_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The main module subroutine.
    ! ----------------------------------------------------------------------
    module subroutine calculate_harmonic_observables_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

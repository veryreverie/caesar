! ======================================================================
! Plots the phonon density of states and dispersion calculated by
!    calculate_harmonic_observables or calculate_anharmonic_observables.
! ======================================================================
module caesar_plot_dos_and_dispersion_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_dos_and_dispersion
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_plot_dos_and_dispersion() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_dos_and_dispersion_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

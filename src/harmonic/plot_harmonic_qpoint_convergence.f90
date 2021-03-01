! ======================================================================
! Plots the convergence of the harmonic free energy calculation
!    w/r/t the q-point grid, as calculated by converge_harmonic_qpoints.
! ======================================================================
module caesar_plot_harmonic_qpoint_convergence_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_harmonic_qpoint_convergence
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_plot_harmonic_qpoint_convergence() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_harmonic_qpoint_convergence_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

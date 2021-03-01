! ======================================================================
! Plots the anharmonic potential along each mode, along with the harmonic
!    and sampled potentials along each mode.
! ======================================================================
module caesar_plot_modes_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_modes
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_plot_modes() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

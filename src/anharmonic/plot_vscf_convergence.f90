! ======================================================================
! Plots the convergence of the VSCF scheme.
! ======================================================================
module caesar_plot_vscf_convergence_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_vscf_convergence
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_plot_vscf_convergence() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_vscf_convergence_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

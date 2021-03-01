! ======================================================================
! Plots the anharmonic potential along multiple modes, along with the harmonic
!    and effective harmonic potentials along those modes mode..
! ======================================================================
module caesar_plot_potential_map_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_potential_map
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_plot_potential_map() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_potential_map_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

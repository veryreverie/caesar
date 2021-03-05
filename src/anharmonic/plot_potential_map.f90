! ======================================================================
! Plots the anharmonic potential along multiple modes, along with the harmonic
!    and effective harmonic potentials along those modes mode..
! ======================================================================
module caesar_plot_potential_map_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: plot_potential_map_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function plot_potential_map_mode() result(output)
      type(ProgramMode) :: output
    end function
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

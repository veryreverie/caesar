! ======================================================================
! Plots the thermodynamic variables calculated by
!    calculate_harmonic_observables or calculate_anharmonic_observables.
! ======================================================================
module caesar_plot_thermodynamic_variables_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: plot_thermodynamic_variables_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function plot_thermodynamic_variables_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_thermodynamic_variables_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

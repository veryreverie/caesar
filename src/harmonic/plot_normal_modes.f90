! ======================================================================
! Plots the results of calculate_normal_modes.
! ======================================================================
module caesar_plot_normal_modes_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: plot_normal_modes_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function plot_normal_modes_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_normal_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

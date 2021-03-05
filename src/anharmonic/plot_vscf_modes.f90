! ======================================================================
! Plots the VSCF potential along each mode.
! ======================================================================
module caesar_plot_vscf_modes_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: plot_vscf_modes_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function plot_vscf_modes_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_vscf_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

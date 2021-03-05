! ======================================================================
! Plots the results of converge_cutoff_and_kpoints.
! ======================================================================
module caesar_plot_cutoff_and_kpoints_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: plot_cutoff_and_kpoints_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function plot_cutoff_and_kpoints_mode() result(output) 
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine plot_cutoff_and_kpoints_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

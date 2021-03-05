! ======================================================================
! Provides routines specific to DFT.
! ======================================================================
! This module is simply an interface for the various DFT modules.
module caesar_dft_module
  use caesar_common_module
  
  use caesar_converge_harmonic_frequencies_module
  use caesar_plot_harmonic_convergence_module
  implicit none
  
  interface
    module function dft_modes() result(output)
      type(ProgramMode), allocatable :: output(:)
    end function
  end interface
end module

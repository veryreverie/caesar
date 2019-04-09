! ======================================================================
! Provides routines specific to DFT.
! ======================================================================
! This module is simply an interface for the various DFT modules.
module dft_module
  use converge_harmonic_frequencies_module
  use plot_harmonic_convergence_module
  implicit none
contains
subroutine startup_dft()
  implicit none
  
  call startup_converge_harmonic_frequencies()
  call startup_plot_harmonic_convergence()
end subroutine
end module

! ======================================================================
! Provides routines specific to CASTEP.
! ======================================================================
! This module is simply an interface for the various castep modules.
module castep_module
  use converge_harmonic_frequencies_module
  use plot_harmonic_convergence_module
  implicit none
contains
subroutine startup_castep()
  implicit none
  
  call startup_converge_harmonic_frequencies()
  call startup_plot_harmonic_convergence()
end subroutine
end module

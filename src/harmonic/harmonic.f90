! ======================================================================
! The operations of Caesar within the harmonic approximation.
! ======================================================================
! This module is simply an interface for the various harmonic modules.
module harmonic_module
  use hartree_to_ev_module
  use snap_to_symmetry_module
  use setup_harmonic_module
  use run_harmonic_module
  use calculate_normal_modes_module
  use plot_normal_modes_module
  use calculate_harmonic_observables_module
  use plot_dos_and_dispersion_module
  use plot_thermodynamic_variables_module
  use converge_qpoint_grid_module
  implicit none
contains
subroutine startup_harmonic()
  implicit none
  
  call startup_hartree_to_ev()
  call startup_snap_to_symmetry()
  call startup_setup_harmonic()
  call startup_run_harmonic()
  call startup_calculate_normal_modes()
  call startup_plot_normal_modes()
  call startup_calculate_harmonic_observables()
  call startup_plot_dos_and_dispersion()
  call startup_plot_thermodynamic_variables()
end subroutine
end module

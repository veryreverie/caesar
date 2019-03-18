! ======================================================================
! The operations of Caesar within the harmonic approximation.
! ======================================================================
! This module is simply an interface for the various harmonic modules.
module harmonic_module
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
end module

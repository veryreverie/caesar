! ======================================================================
! The operations of Caesar beyond the harmonic approximation.
! ======================================================================
! This module is simply an interface for the various harmonic modules.
module anharmonic_module
  use setup_anharmonic_module
  use run_anharmonic_module
  use calculate_potential_module
  use map_anharmonic_modes_module
  use plot_anharmonic_modes_module
  use map_potential_module
  use plot_potential_map_module
  use map_vscf_modes_module
  use plot_vscf_modes_module
  use calculate_anharmonic_observables_module
  use plot_vscf_convergence_module
  use plot_vscf_states_module
  use states_module
  use potentials_module
  implicit none
contains
subroutine startup_anharmonic()
  implicit none
  
  call startup_states()
  call startup_potentials()
  call startup_setup_anharmonic()
  call startup_run_anharmonic()
  call startup_calculate_potential()
  call startup_map_anharmonic_modes()
  call startup_plot_anharmonic_modes()
  call startup_map_potential()
  call startup_plot_potential_map()
  call startup_map_vscf_modes()
  call startup_plot_vscf_modes()
  call startup_calculate_anharmonic_observables()
  call startup_plot_vscf_convergence()
  call startup_plot_vscf_states()
end subroutine
end module

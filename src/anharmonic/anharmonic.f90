! ======================================================================
! The operations of Caesar beyond the harmonic approximation.
! ======================================================================
! This module is simply an interface for the various harmonic modules.
module caesar_anharmonic_module
  use caesar_setup_anharmonic_module
  use caesar_run_anharmonic_module
  use caesar_calculate_potential_module
  use caesar_map_modes_module
  use caesar_plot_modes_module
  use caesar_map_potential_module
  use caesar_plot_potential_map_module
  use caesar_map_vscf_modes_module
  use caesar_plot_vscf_modes_module
  use caesar_calculate_anharmonic_observables_module
  use caesar_plot_vscf_convergence_module
  use caesar_plot_vscf_states_module
  use caesar_states_module
  use caesar_potentials_module
  implicit none
contains
subroutine startup_anharmonic()
  implicit none
  
  call startup_states()
  call startup_potentials()
  call startup_setup_anharmonic()
  call startup_run_anharmonic()
  call startup_calculate_potential()
  call startup_map_modes()
  call startup_plot_modes()
  call startup_map_potential()
  call startup_plot_potential_map()
  call startup_map_vscf_modes()
  call startup_plot_vscf_modes()
  call startup_calculate_anharmonic_observables()
  call startup_plot_vscf_convergence()
  call startup_plot_vscf_states()
end subroutine
end module

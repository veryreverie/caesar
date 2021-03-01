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
  
  interface
    module subroutine startup_anharmonic()
    end subroutine
  end interface
end module

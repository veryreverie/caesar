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
  use plot_vscf_states_module
  use states_module
  use potentials_module
  implicit none
contains
subroutine startup_anharmonic()
  implicit none
  
  call startup_states()
  call startup_potentials()
end subroutine
end module

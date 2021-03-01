submodule (caesar_anharmonic_module) caesar_anharmonic_submodule
contains
module procedure startup_anharmonic
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
end procedure
end submodule

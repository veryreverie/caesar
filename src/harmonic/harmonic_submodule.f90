submodule (caesar_harmonic_module) caesar_harmonic_submodule
contains
module procedure startup_harmonic
  call startup_hartree_to_ev()
  call startup_snap_to_symmetry()
  call startup_setup_harmonic()
  call startup_run_harmonic()
  call startup_calculate_normal_modes()
  call startup_read_normal_modes()
  call startup_plot_normal_modes()
  call startup_calculate_harmonic_observables()
  call startup_plot_dos_and_dispersion()
  call startup_plot_thermodynamic_variables()
  call startup_converge_harmonic_qpoints()
  call startup_plot_harmonic_qpoint_convergence()
end procedure
end submodule

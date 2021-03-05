submodule (caesar_harmonic_module) caesar_harmonic_submodule
contains
module procedure harmonic_modes
  output = [                                   &
     & hartree_to_ev_mode(),                   &
     & snap_to_symmetry_mode(),                &
     & setup_harmonic_mode(),                  &
     & run_harmonic_mode(),                    &
     & calculate_normal_modes_mode(),          &
     & read_normal_modes_mode(),               &
     & plot_normal_modes_mode(),               &
     & calculate_harmonic_observables_mode(),  &
     & plot_dos_and_dispersion_mode(),         &
     & plot_thermodynamic_variables_mode(),    &
     & converge_harmonic_qpoints_mode(),       &
     & plot_harmonic_qpoint_convergence_mode() ]
end procedure
end submodule

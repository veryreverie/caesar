submodule (caesar_anharmonic_module) caesar_anharmonic_submodule
contains

module procedure anharmonic_modes
  output = [                                    &
     & setup_anharmonic_mode(),                 &
     & run_anharmonic_mode(),                   &
     & calculate_potential_mode(),              &
     & map_modes_mode(),                        &
     & plot_modes_mode(),                       &
     & map_potential_mode(),                    &
     & plot_potential_map_mode(),               &
     & map_vscf_modes_mode(),                   &
     & plot_vscf_modes_mode(),                  &
     & calculate_anharmonic_observables_mode(), &
     & plot_vscf_convergence_mode(),            &
     & plot_vscf_states_mode()                  ]
end procedure
end submodule

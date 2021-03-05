submodule (caesar_dft_module) caesar_dft_submodule
contains
module procedure dft_modes
  output = [                                 &
     & converge_harmonic_frequencies_mode(), &
     & plot_harmonic_convergence_mode()      ]
end procedure
end submodule

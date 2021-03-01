submodule (caesar_dft_module) caesar_dft_submodule
contains
module procedure startup_dft
  call startup_converge_harmonic_frequencies()
  call startup_plot_harmonic_convergence()
end procedure
end submodule

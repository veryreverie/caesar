submodule (caesar_polynomial_module) caesar_polynomial_submodule
contains

module procedure startup_polynomial
  call startup_polynomial_stress()
  call startup_polynomial_potential()
end procedure
end submodule

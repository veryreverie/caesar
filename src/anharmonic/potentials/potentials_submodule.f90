submodule (caesar_potentials_module) caesar_potentials_submodule
contains
module procedure startup_potentials
  call startup_polynomial()
  call startup_potential_example()
end procedure
end submodule

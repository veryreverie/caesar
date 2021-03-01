submodule (caesar_basic_atom_module) caesar_basic_atom_submodule
  use caesar_atom_module
contains

module procedure new_BasicAtom
  output%species            = species
  output%mass               = mass
  output%cartesian_position = cartesian_position
end procedure
end submodule

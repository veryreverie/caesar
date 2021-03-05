submodule (caesar_basis_states_module) caesar_basis_states_types_submodule
  use caesar_effective_harmonic_module
  use caesar_wavevector_states_module
contains

module procedure types_BasisStates
  type(HarmonicStates)   :: type_1
  type(WavevectorStates) :: type_2
  
  output = [ BasisStatesPointer(type_1), &
           & BasisStatesPointer(type_2)  ]
end procedure
end submodule

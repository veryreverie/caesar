submodule (caesar_subspace_state_module) caesar_subspace_state_types_submodule
  !use caesar_monomial_state_real_module
  !use caesar_monomial_state_complex_module
  use caesar_harmonic_state_real_module
  use caesar_harmonic_state_complex_module
contains

module procedure types_SubspaceState
  type(HarmonicStateReal)    :: type_1
  type(HarmonicStateComplex) :: type_2
  
  output = [ SubspaceStatePointer(type_1), &
           & SubspaceStatePointer(type_2)  ]
end procedure
end submodule

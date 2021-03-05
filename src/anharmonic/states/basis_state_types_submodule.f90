submodule (caesar_basis_state_module) caesar_basis_state_types_submodule
  use caesar_wavevector_state_module
contains

module procedure types_BasisState
  type(WavevectorState) :: type_1
  
  output = [BasisStatePointer(type_1)]
end procedure
end submodule

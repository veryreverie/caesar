submodule (caesar_abstract_classes_module) caesar_abstract_state_types_submodule
  use caesar_effective_harmonic_module
  
  use caesar_wavevector_basis_module
  use caesar_full_subspace_basis_module
  use caesar_split_qpoints_basis_module
contains

module procedure types_SubspaceBasis
  type(HarmonicBasis)     :: type_1
  type(WavevectorBasis)   :: type_2
  type(FullSubspaceBasis) :: type_3
  type(SplitQpointsBasis) :: type_4
  
  output = [ SubspaceBasisPointer(type_1), &
           & SubspaceBasisPointer(type_2), &
           & SubspaceBasisPointer(type_3), &
           & SubspaceBasisPointer(type_4)  ]
end procedure
end submodule

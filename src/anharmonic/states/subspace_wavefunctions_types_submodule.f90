submodule (caesar_subspace_wavefunctions_module) caesar_subspace_wavefunctions_types_submodule
  use caesar_full_subspace_wavefunctions_module
  use caesar_split_qpoints_wavefunctions_module
contains

module procedure types_SubspaceWavefunctions
  type(FullSubspaceWavefunctions) :: type_1
  type(SplitQpointsWavefunctions) :: type_2
  
  output = [ SubspaceWavefunctionsPointer(type_1), &
           & SubspaceWavefunctionsPointer(type_2)  ]
end procedure
end submodule

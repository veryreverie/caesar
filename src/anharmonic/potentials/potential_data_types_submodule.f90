submodule (caesar_potential_data_module) caesar_potential_data_types_submodule
  use caesar_polynomial_module
  
  use caesar_potential_example_module
contains

module procedure types_PotentialData
  type(PolynomialPotential)  :: type_1
  type(PotentialDataExample) :: type_2
  
  output = [ PotentialPointer(type_1), &
           & PotentialPointer(type_2)  ]
end procedure
end submodule

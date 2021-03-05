submodule (caesar_stress_data_module) caesar_stress_data_types_submodule
  use caesar_polynomial_module
contains

module procedure types_StressData
  type(PolynomialStress) :: type_1
  
  output = [StressPointer(type_1)]
end procedure
end submodule

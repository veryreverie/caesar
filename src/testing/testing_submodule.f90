submodule (caesar_testing_module) caesar_testing_submodule
contains

module procedure testing_modes
  output = [                          &
     & check_counter_mode(),          &
     & update_basis_functions_mode(), &
     & test_mode()                    ]
end procedure
end submodule

submodule (caesar_testing_module) caesar_testing_submodule
contains

module procedure startup_testing
  call startup_check_counter()
  call startup_update_basis_functions()
  call startup_test()
end procedure
end submodule

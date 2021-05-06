cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(testing)

add_library(testing_modules atom_mapping.f90 displacement_patterns.f90 linear_algebra_test.f90 structure_test.f90 check_counter.f90 update_basis_functions.f90 test.f90 testing.f90)
add_library(testing_submodules atom_mapping_submodule.f90 displacement_patterns_submodule.f90 linear_algebra_test_submodule.f90 structure_test_submodule.f90 check_counter_submodule.f90 update_basis_functions_submodule.f90 test_submodule.f90 testing_submodule.f90)
target_link_libraries(testing_submodules testing_modules)
add_dependencies(testing_submodules testing_modules)

# Link against required modules.
foreach(dir IN ITEMS common harmonic anharmonic anharmonic_common)
  target_link_libraries(testing_modules ${dir}_modules)
  target_link_libraries(testing_submodules ${dir}_submodules)
  add_dependencies(testing_modules ${dir}_modules)
endforeach()

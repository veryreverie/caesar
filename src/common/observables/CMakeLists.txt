cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(observables)

# Normal mode modules.
add_library(observables_modules thermodynamic_data.f90 qpoint_path.f90 harmonic_properties.f90 observables.f90)
add_library(observables_submodules thermodynamic_data_submodule.f90 qpoint_path_submodule.f90 harmonic_properties_submodule.f90)
target_link_libraries(observables_submodules observables_modules)

# Link against required module.
foreach(dir IN ITEMS utils structure normal_mode dynamical_matrices harmonic_stress)
  target_link_libraries(observables_modules ${dir}_modules)
  target_link_libraries(observables_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs thermodynamic_data_test.pf qpoint_path_test.pf harmonic_properties_test.pf)
  add_pfunit_ctest(observables_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(observables_tests observables_submodules)
endif()

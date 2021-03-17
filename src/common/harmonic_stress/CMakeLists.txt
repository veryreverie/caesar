cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(harmonic_stress)

# Normal mode modules.
add_library(harmonic_stress_modules stress_hessian.f90 stress_dynamical_matrix.f90 harmonic_stress_conversion.f90 harmonic_stress.f90)
add_library(harmonic_stress_submodules stress_hessian_submodule.f90 stress_dynamical_matrix_submodule.f90 harmonic_stress_conversion_submodule.f90)
target_link_libraries(harmonic_stress_submodules harmonic_stress_modules)
add_dependencies(harmonic_stress_submodules harmonic_stress_modules)

# Link against required module.
foreach(dir IN ITEMS utils structure normal_mode dynamical_matrices)
  target_link_libraries(harmonic_stress_modules ${dir}_modules)
  target_link_libraries(harmonic_stress_submodules ${dir}_submodules)
  add_dependencies(harmonic_stress_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs stress_hessian_test.pf stress_dynamical_matrix_test.pf harmonic_stress_conversion_test.pf)
  add_pfunit_ctest(harmonic_stress_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(harmonic_stress_tests harmonic_stress_submodules)
endif()

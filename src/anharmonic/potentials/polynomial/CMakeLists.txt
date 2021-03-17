cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(polynomial)

add_library(polynomial_modules polynomial_dynamical_matrices.f90 polynomial_interpolator.f90 polynomial_symmetry.f90 basis_function.f90 sampling_points.f90 coupling_basis_functions.f90 stress_basis_function.f90 coupling_stress_basis_functions.f90 vscf_rvector.f90 vscf_rvectors.f90 sample_result.f90 sample_results.f90 polynomial_stress.f90 polynomial_potential.f90 permutation.f90 polynomial.f90)
add_library(polynomial_submodules polynomial_dynamical_matrices_submodule.f90 polynomial_interpolator_submodule.f90 polynomial_symmetry_submodule.f90 basis_function_submodule.f90 sampling_points_submodule.f90 coupling_basis_functions_submodule.f90 stress_basis_function_submodule.f90 coupling_stress_basis_functions_submodule.f90 vscf_rvector_submodule.f90 vscf_rvectors_submodule.f90 sample_result_submodule.f90 polynomial_stress_submodule.f90 polynomial_potential_submodule.f90 permutation_submodule.f90)
target_link_libraries(polynomial_submodules polynomial_modules)
add_dependencies(polynomial_submodules polynomial_modules)

# Link against required modules.
foreach(dir IN ITEMS common harmonic anharmonic_common states)
  target_link_libraries(polynomial_modules ${dir}_modules)
  target_link_libraries(polynomial_submodules ${dir}_submodules)
  add_dependencies(polynomial_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs polynomial_dynamical_matrices_test.pf polynomial_interpolator_test.pf polynomial_symmetry_test.pf basis_function_test.pf sampling_points_test.pf coupling_basis_functions_test.pf stress_basis_function_test.pf coupling_stress_basis_functions_test.pf vscf_rvector_test.pf vscf_rvectors_test.pf sample_result_test.pf polynomial_stress_test.pf polynomial_potential_test.pf permutation_test.pf)
  add_pfunit_ctest(polynomial_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(polynomial_tests polynomial_submodules)
endif()

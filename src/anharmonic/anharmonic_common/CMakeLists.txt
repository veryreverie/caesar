cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(anharmonic_common)

add_library(anharmonic_common_modules interpolated_supercell.f90 max_displacement.f90 subspace_coupling.f90 subspace_monomial.f90 degenerate_symmetry.f90 anharmonic_data.f90 subspace_wavefunctions.f90 stress_prefactors.f90 sparse_monomial.f90 subspace_state.f90 subspace_braket.f90 basis_state.f90 expectation_cache.f90 basis_states.f90 pulay.f90 abstract_classes.f90 braket.f90 stress_data.f90 potential_data.f90 anharmonic_common.f90)
add_library(anharmonic_common_submodules interpolated_supercell_submodule.f90 max_displacement_submodule.f90 subspace_coupling_submodule.f90 subspace_monomial_submodule.f90 degenerate_symmetry_submodule.f90 anharmonic_data_submodule.f90 subspace_wavefunctions_submodule.f90 stress_prefactors_submodule.f90 sparse_monomial_submodule.f90 subspace_state_submodule.f90 basis_state_submodule.f90 expectation_cache_submodule.f90 basis_states_submodule.f90 pulay_submodule.f90 abstract_classes_submodule.f90 braket_submodule.f90 stress_data_submodule.f90 potential_data_submodule.f90)
target_link_libraries(anharmonic_common_submodules anharmonic_common_modules)
target_include_directories(anharmonic_common_modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
foreach(dir IN ITEMS integration)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(anharmonic_common_modules ${dir}_modules)
  target_link_libraries(anharmonic_common_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS common harmonic)
  target_link_libraries(anharmonic_common_modules ${dir}_modules)
  target_link_libraries(anharmonic_common_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs interpolated_supercell_test.pf max_displacement_test.pf subspace_coupling_test.pf subspace_monomial_test.pf degenerate_symmetry_test.pf anharmonic_data_test.pf subspace_wavefunctions_test.pf stress_prefactors_test.pf sparse_monomial_test.pf subspace_state_test.pf subspace_braket_test.pf basis_state_test.pf expectation_cache_test.pf basis_states_test.pf pulay_test.pf abstract_classes_test.pf braket_test.pf stress_data_test.pf potential_data_test.pf)
  add_pfunit_ctest(anharmonic_common_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(anharmonic_common_tests anharmonic_common_submodules)
endif()

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(states)

add_library(states_modules harmonic_state_1d.f90 harmonic_state_2d.f90 harmonic_state_real.f90 harmonic_state_complex.f90 harmonic_braket_real.f90 harmonic_braket_complex.f90 coupled_states.f90 density_matrix.f90 wavevector_state.f90 wavevector_states.f90 wavevector_basis.f90 full_subspace_wavefunctions.f90 split_qpoints_wavefunctions.f90 core_shell_thermodynamics.f90 calculate_weights.f90 full_subspace_basis.f90 split_qpoints_basis.f90 states.f90)
add_library(states_submodules harmonic_state_1d_submodule.f90 harmonic_state_2d_submodule.f90 harmonic_state_real_submodule.f90 harmonic_state_complex_submodule.f90 harmonic_braket_real_submodule.f90 harmonic_braket_complex_submodule.f90 coupled_states_submodule.f90 density_matrix_submodule.f90 wavevector_state_submodule.f90 wavevector_states_submodule.f90 wavevector_basis_submodule.f90 full_subspace_wavefunctions_submodule.f90 split_qpoints_wavefunctions_submodule.f90 core_shell_thermodynamics_submodule.f90 calculate_weights_submodule.f90 full_subspace_basis_submodule.f90 split_qpoints_basis_submodule.f90 subspace_state_types_submodule.f90 basis_state_types_submodule.f90 basis_states_types_submodule.f90 abstract_states_types_submodule.f90 subspace_wavefunctions_types_submodule.f90)
target_link_libraries(states_submodules states_modules)

# Include subsidiary modules.
foreach(dir IN ITEMS effective_harmonic)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(states_modules ${dir}_modules)
  target_link_libraries(states_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS common harmonic anharmonic_common)
  target_link_libraries(states_modules ${dir}_modules)
  target_link_libraries(states_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs harmonic_state_1d_test.pf harmonic_state_2d_test.pf harmonic_state_real_test.pf harmonic_state_complex_test.pf harmonic_braket_real_test.pf harmonic_braket_complex_test.pf coupled_states_test.pf density_matrix_test.pf wavevector_state_test.pf wavevector_states_test.pf wavevector_basis_test.pf full_subspace_wavefunctions_test.pf split_qpoints_wavefunctions_test.pf core_shell_thermodynamics_test.pf calculate_weights_test.pf full_subspace_basis_test.pf split_qpoints_basis_test.pf)
  add_pfunit_ctest(states_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(states_tests states_submodules)
endif()

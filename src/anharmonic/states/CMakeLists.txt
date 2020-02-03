cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(states)

add_library(states monomial_state_1d.f90 monomial_state_2d.f90 harmonic_state_1d.f90 harmonic_state_2d.f90 monomial_state_real.f90 monomial_state_complex.f90 harmonic_state_real.f90 harmonic_state_complex.f90 coupled_states.f90 density_matrix.f90 wavevector_state.f90 wavevector_states.f90 wavevector_basis.f90 full_subspace_wavefunctions.f90 split_qpoints_wavefunctions.f90 core_shell_thermodynamics.f90 calculate_weights.f90 full_subspace_basis.f90 split_qpoints_basis.f90 states.f90)
target_include_directories(states PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/effective_harmonic effective_harmonic/)
target_link_libraries(states effective_harmonic)

# Link against other modules.
target_link_libraries(states harmonic)
target_link_libraries(states utils)
target_link_libraries(states common)
target_link_libraries(states anharmonic_common)

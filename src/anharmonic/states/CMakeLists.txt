cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(states)

add_library(states state_helper.f90 monomial_state.f90 polynomial_state.f90 harmonic_state.f90 braket.f90 state_conversion.f90 coupled_states.f90 wavevector_basis.f90 full_subspace_basis_and_states.f90 subspace_basis_pointer.f90 subspace_states_pointer.f90 states.f90)
target_include_directories(states PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(states harmonic)
target_link_libraries(states utils)
target_link_libraries(states common)
target_link_libraries(states anharmonic_common)

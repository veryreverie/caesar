cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic_common)

add_library(anharmonic_common degeneracy.f90 subspace_coupling.f90 subspace_monomial.f90 mode_coupling.f90 mode_monomial.f90 degenerate_symmetry.f90 anharmonic_data.f90 single_mode_state.f90 single_harmonic_basis.f90 single_harmonic_state.f90 potential.f90 potential_pointer.f90 effective_frequency.f90 anharmonic_common.f90)
target_include_directories(anharmonic_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(anharmonic_common harmonic)
target_link_libraries(anharmonic_common utils)
target_link_libraries(anharmonic_common common)

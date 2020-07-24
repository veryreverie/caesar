cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic_common)

add_library(anharmonic_common interpolated_supercell.f90 max_displacement.f90 subspace_coupling.f90 subspace_monomial.f90 degenerate_symmetry.f90 anharmonic_data.f90 subspace_wavefunctions.f90 stress_prefactors.f90 sparse_monomial.f90 subspace_state.f90 subspace_braket.f90 basis_state.f90 expectation_cache.f90 basis_states.f90 pulay.f90 abstract_classes.f90 braket.f90 stress_data.f90 potential_data.f90 anharmonic_common.f90)
target_include_directories(anharmonic_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/integration integration/)
target_link_libraries(anharmonic_common integration)

# Link against other modules.
target_link_libraries(anharmonic_common harmonic)
target_link_libraries(anharmonic_common utils)
target_link_libraries(anharmonic_common common)

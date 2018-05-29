cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(shared)

add_library(shared degeneracy.f90 coupled_subspaces.f90 subspace_monomial.f90 coupled_modes.f90 mode_monomial.f90 degenerate_symmetry.f90 single_mode_state.f90 subspace_product_state.f90 potential.f90 shared.f90)
target_include_directories(shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(shared harmonic)
target_link_libraries(shared utils)
target_link_libraries(shared common)

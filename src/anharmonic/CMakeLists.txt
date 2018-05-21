cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic)

add_library(anharmonic degeneracy.f90 coupled_subspaces.f90 subspace_monomial.f90 coupled_modes.f90 mode_monomial.f90 degenerate_symmetry.f90 basis_function.f90 basis_functions.f90 sampling_points.f90 vscf_rvectors.f90 setup_anharmonic.f90)
target_include_directories(anharmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/harmonic harmonic/)
target_link_libraries(anharmonic harmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(anharmonic utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/common common/)
target_link_libraries(anharmonic common)

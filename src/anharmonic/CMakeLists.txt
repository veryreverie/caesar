cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic)

add_library(anharmonic generate_subspace_potentials.f90 mode_map.f90 vscf.f90 interpolation.f90 setup_anharmonic.f90 run_anharmonic.f90 calculate_potential.f90 map_anharmonic_modes.f90 plot_anharmonic_modes.f90 map_potential.f90 plot_potential_map.f90 map_vscf_modes.f90 plot_vscf_modes.f90 calculate_anharmonic_observables.f90 plot_vscf_states.f90 plot_vscf_convergence.f90 anharmonic.f90)
target_include_directories(anharmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/anharmonic_common anharmonic_common/)
target_link_libraries(anharmonic anharmonic_common)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/states states/)
target_link_libraries(anharmonic states)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/potentials potentials/)
target_link_libraries(anharmonic potentials)

# Link against other modules.
target_link_libraries(anharmonic harmonic)
target_link_libraries(anharmonic utils)
target_link_libraries(anharmonic common)

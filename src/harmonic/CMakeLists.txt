cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(harmonic)

add_library(harmonic generate_supercells.f90 hartree_to_ev.f90 snap_to_symmetry.f90 construct_supercell_hessian.f90 harmonic_data.f90 setup_harmonic.f90 run_harmonic.f90 calculate_normal_modes.f90 read_normal_modes.f90 calculate_harmonic_observables.f90 plot_normal_modes.f90 plot_dos_and_dispersion.f90 plot_thermodynamic_variables.f90 converge_qpoint_grid.f90 converge_harmonic_qpoints.f90 harmonic.f90)
target_include_directories(harmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(harmonic utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/common common/)
target_link_libraries(harmonic common)

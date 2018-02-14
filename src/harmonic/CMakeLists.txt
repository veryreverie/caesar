cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(harmonic)

add_library(harmonic unique_directions.f90 min_images.f90 generate_supercells.f90 dynamical_matrix.f90 lte.f90 force_constants.f90 hartree_to_ev.f90 setup_harmonic.f90 run_harmonic.f90 calculate_normal_modes.f90 calculate_dos_and_dispersion.f90 plot_normal_modes.f90 plot_dos_and_dispersion.f90)
target_include_directories(harmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(harmonic utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/common common/)
target_link_libraries(harmonic common)

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(castep)

add_library(castep castep_output_file.f90 converge_cutoff_and_kpoints.f90 plot_cutoff_and_kpoints.f90 converge_harmonic_frequencies.f90 plot_harmonic_convergence.f90)
target_include_directories(castep PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(castep common)

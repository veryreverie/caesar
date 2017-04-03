cmake_minimum_required (version 2.8 fatal_error)
include(../utils/cmake_include_guard)
include_guard(harmonic)

add_library(harmonic calculate_symmetry_group.f90 unique_directions.f90 construct_supercell.f90 min_images.f90 generate_supercells.f90 lte.f90 hartree_to_eV.f90 setup_harmonic.f90 run_harmonic.f90 lte_harmonic.f90)
target_include_directories(harmonic public ${cmake_current_source_dir})

add_subdirectory(utils utils/)
target_link_libraries(harmonic utils)

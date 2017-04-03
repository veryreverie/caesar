cmake_minimum_required (version 2.8 fatal_error)
include(../cmake_include_guard)
include_guard(utils)

add_library(utils constants.f90 string.f90 file.f90 utils.f90 linear_algebra.f90 algebra.f90 group.f90 supercell.f90 kpoints.f90 moller_plesset.f90 structure.f90 dft_input_file.f90 dft_output_file.f90 structure_to_dft.f90 calculate_symmetry.f90 bands.f90 displacement_patterns.f90)
target_include_directories(utils public ${cmake_current_source_dir})

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(castep)

add_library(castep castep_wrapper.f90 phonon_file.f90 castep.f90)

# Link against required modules.
target_link_libraries(castep utils)
target_link_libraries(castep structure)
target_link_libraries(castep normal_mode)
target_link_libraries(castep dynamical_matrices)
target_link_libraries(castep electronic_structure_data)
target_link_libraries(castep electronic_structure_common)

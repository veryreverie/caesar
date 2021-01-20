cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(vasp)

add_library(vasp vasp_wrapper.f90 vasp.f90)

# Link against required modules.
target_link_libraries(vasp utils)
target_link_libraries(vasp structure)
target_link_libraries(vasp normal_mode)
target_link_libraries(vasp dynamical_matrices)
target_link_libraries(vasp electronic_structure_data)
target_link_libraries(vasp electronic_structure_common)

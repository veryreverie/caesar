cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(qe)

add_library(qe qe_fc.f90 qe_pw.f90 qe.f90)

# Link against required modules.
target_link_libraries(qe utils)
target_link_libraries(qe structure)
target_link_libraries(qe normal_mode)
target_link_libraries(qe dynamical_matrices)
target_link_libraries(qe electronic_structure_data)
target_link_libraries(qe electronic_structure_common)

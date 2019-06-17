cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(dynamical_matrices)

# Normal mode modules.
add_library(dynamical_matrices unique_directions.f90 min_images.f90 construct_hessian.f90 dynamical_matrix.f90 dynamical_matrices.f90)

# Link against required module.
target_link_libraries(dynamical_matrices utils)
target_link_libraries(dynamical_matrices structure)
target_link_libraries(dynamical_matrices normal_mode)
target_link_libraries(dynamical_matrices electronic_structure)

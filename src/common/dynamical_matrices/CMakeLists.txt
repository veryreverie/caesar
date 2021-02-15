cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(dynamical_matrices)

# Normal mode modules.
add_library(dynamical_matrices unique_directions.f90 min_images.f90 dynamical_matrix.f90 calculate_dynamical_matrices.f90 dynamical_matrices.f90)

# Link against required module.
target_link_libraries(dynamical_matrices utils)
target_link_libraries(dynamical_matrices structure)
target_link_libraries(dynamical_matrices normal_mode)

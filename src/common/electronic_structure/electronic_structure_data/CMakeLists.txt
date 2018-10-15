cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(electronic_structure_data)

add_library(electronic_structure_data linear_response.f90 electronic_structure_data.f90)

# Link against required modules.
target_link_libraries(electronic_structure_data utils)
target_link_libraries(electronic_structure_data structure)
target_link_libraries(electronic_structure_data normal_mode)

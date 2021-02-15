cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(xyz)

add_library(xyz xyz.f90)

# Link against required modules.
target_link_libraries(xyz utils)
target_link_libraries(xyz structure)

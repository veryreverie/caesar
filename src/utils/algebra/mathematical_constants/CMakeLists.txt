cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(mathematical_constants)

add_library(mathematical_constants mathematical_constants.f90)

# Link against required modules.
target_link_libraries(mathematical_constants precision)

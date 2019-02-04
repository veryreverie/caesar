cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(clock)

add_library(clock clock.f90)

# Link against required modules.
target_link_libraries(clock precision)

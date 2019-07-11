cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(profiling)

add_library(profiling clock.f90 getpid.c memory_checker.f90 profiler.f90 profiling.f90)

# Link against required modules.
target_link_libraries(profiling precision)

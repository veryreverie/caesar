cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(harmonic_stress)

# Normal mode modules.
add_library(harmonic_stress harmonic_stress.f90)

# Link against required module.
target_link_libraries(harmonic_stress utils)
target_link_libraries(harmonic_stress structure)

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(observables)

# Normal mode modules.
add_library(observables harmonic_thermodynamics.f90 observables.f90)

# Link against required module.
target_link_libraries(observables utils)
target_link_libraries(observables structure)
target_link_libraries(observables normal_mode)
target_link_libraries(observables electronic_structure)

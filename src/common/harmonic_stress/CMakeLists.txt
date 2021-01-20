cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(harmonic_stress)

# Normal mode modules.
add_library(harmonic_stress stress_hessian.f90 stress_dynamical_matrix.f90 harmonic_stress_conversion.f90 harmonic_stress.f90)

# Link against required module.
target_link_libraries(harmonic_stress utils)
target_link_libraries(harmonic_stress structure)
target_link_libraries(harmonic_stress normal_mode)
target_link_libraries(harmonic_stress dynamical_matrices)

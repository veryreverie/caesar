cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(normal_mode)

# Normal mode modules.
add_library(normal_mode cartesian_vector.f90 cartesian_displacement.f90 cartesian_force.f90 real_mode.f90 complex_mode.f90 real_single_mode_vector.f90 complex_single_mode_vector.f90 real_mode_vector.f90 complex_mode_vector.f90 real_mode_displacement.f90 complex_mode_displacement.f90 real_mode_force.f90
complex_mode_force.f90 complex_polynomial.f90 real_polynomial.f90 real_complex_conversion.f90 complex_mode_symmetry.f90 normal_mode.f90)

# Link against required module.
target_link_libraries(normal_mode utils)
target_link_libraries(normal_mode structure)

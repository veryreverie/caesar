cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(normal_mode)

# Normal mode modules.
add_library(normal_mode cartesian_displacement.f90 cartesian_force.f90 cartesian_hessian.f90 mass_weighted_displacement.f90 mass_weighted_force.f90 real_mode.f90 complex_mode.f90 real_single_mode_displacement.f90 real_single_mode_force.f90 complex_single_mode_displacement.f90 complex_single_mode_force.f90 real_mode_displacement.f90 real_mode_force.f90 complex_mode_displacement.f90 complex_mode_force.f90 complex_polynomial.f90 real_polynomial.f90 paired_polynomial.f90 real_complex_conversion.f90 complex_mode_symmetry.f90 degenerate_subspace.f90 process_modes.f90 normal_mode.f90)

# Link against required module.
target_link_libraries(normal_mode utils)
target_link_libraries(normal_mode structure)

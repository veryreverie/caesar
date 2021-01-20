cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(observables)

# Normal mode modules.
add_library(observables thermodynamic_data.f90 qpoint_path.f90 harmonic_properties.f90 observables.f90)

# Link against required module.
target_link_libraries(observables utils)
target_link_libraries(observables structure)
target_link_libraries(observables normal_mode)
target_link_libraries(observables dynamical_matrices)
target_link_libraries(observables harmonic_stress)

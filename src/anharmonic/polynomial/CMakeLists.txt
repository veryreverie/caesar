cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(polynomial)

add_library(polynomial basis_function.f90 basis_functions.f90 sampling_points.f90 vscf_rvector.f90 vscf_rvectors.f90 polynomial_potential.f90 polynomial.f90)
target_include_directories(polynomial PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(polynomial harmonic)
target_link_libraries(polynomial utils)
target_link_libraries(polynomial common)
target_link_libraries(polynomial anharmonic_common)

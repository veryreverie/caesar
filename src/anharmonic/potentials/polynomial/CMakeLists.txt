cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(polynomial)

add_library(polynomial polynomial_interpolator.f90 basis_function.f90 sampling_points.f90 coupling_basis_functions.f90 stress_basis_function.f90 coupling_stress_basis_functions.f90 vscf_rvector.f90 vscf_rvectors.f90 sample_result.f90 sample_results.f90 fit_coefficients.f90 fit_stress_coefficients.f90 polynomial_interpolation.f90 polynomial_stress.f90 polynomial_potential.f90 permutation.f90 polynomial.f90)
target_include_directories(polynomial PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(polynomial harmonic)
target_link_libraries(polynomial utils)
target_link_libraries(polynomial common)
target_link_libraries(polynomial states)
target_link_libraries(polynomial anharmonic_common)

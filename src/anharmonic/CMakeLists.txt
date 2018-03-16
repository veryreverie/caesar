cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic)

#add_library(anharmonic mode_vector.f90 sampling_points.f90 coupling.f90 single_mode_states.f90 harmonic_states.f90 vscf_states.f90 product_states.f90 monomial.f90 potential.f90 grid_types.f90 scf.f90 perturbation.f90 setup_anharmonic.f90 run_anharmonic.f90 calculate_anharmonic.f90)
add_library(anharmonic single_mode_displacement.f90 mode_displacement.f90 univariate.f90 monomial.f90 polynomial.f90 degeneracy.f90 coupling.f90 degenerate_symmetry.f90 setup_anharmonic.f90 basis_function.f90)# real_mode.f90)# displacements_and_forces.f90 grid_types.f90 sampling_points.f90 single_mode_states.f90 harmonic_states.f90 vscf_states.f90 product_states.f90 potential.f90 scf.f90 perturbation.f90 run_anharmonic.f90 calculate_anharmonic.f90)
target_include_directories(anharmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/harmonic harmonic/)
target_link_libraries(anharmonic harmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(anharmonic utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/common common/)
target_link_libraries(anharmonic common)

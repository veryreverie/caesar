cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(quadratic)

add_library(quadratic calculate_anharmonic_correction.f90 calculate_gap.f90 quadratic_spline.f90 vscf_1d.f90 setup_quadratic.f90 run_quadratic.f90 anharmonic.f90 bs_quadratic.f90)
target_include_directories(harmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(quadratic utils)

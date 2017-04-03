cmake_minimum_required (version 2.8 fatal_error)
include(../utils/cmake_include_guard)
include_guard(quadratic)

add_library(quadratic mapping.f90 calculate_anharmonic.f90 calculate_gap.f90 quadratic_spline.f90 vscf_1d.f90 setup_quadratic.f90 run_quadratic.f90 anharmonic.f90 bs_quadratic.f90)
target_include_directories(harmonic public ${cmake_current_source_dir})

add_subdirectory(utils utils/)
target_link_libraries(quadratic utils)

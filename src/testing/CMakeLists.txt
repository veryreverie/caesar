cmake_minimum_required (version 2.8 fatal_error)
include(../utils/cmake_include_guard)
include_guard(testing)

add_library(testing atom_mapping.f90 test_copy_harmonic.f90 test_lte.f90 test_copy_quadratic.f90)
target_include_directories(testing public ${cmake_current_source_dir})

add_subdirectory(utils utils/)
target_link_libraries(testing utils)

add_subdirectory(harmonic harmonic/)
target_link_libraries(testing harmonic)

add_subdirectory(quadratic quadratic/)
target_link_libraries(testing quadratic)

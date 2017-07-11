cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(testing)

add_library(testing atom_mapping.f90 displacement_patterns.f90 linear_algebra_test.f90 structure_test.f90 setup_harmonic_test.f90 lte_harmonic_test.f90 test_copy_quadratic.f90 test.f90)
target_include_directories(testing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(testing utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/harmonic harmonic/)
target_link_libraries(testing harmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/quadratic quadratic/)
target_link_libraries(testing quadratic)

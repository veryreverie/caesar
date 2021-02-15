cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(effective_harmonic)

add_library(effective_harmonic harmonic_states.f90 harmonic_basis.f90 effective_harmonic.f90)
target_include_directories(effective_harmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(effective_harmonic utils)
target_link_libraries(effective_harmonic common)
target_link_libraries(effective_harmonic anharmonic_common)

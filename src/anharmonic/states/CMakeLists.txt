cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(states)

add_library(states subspace_state.f90 subspace_basis.f90 states.f90)
target_include_directories(states PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
target_link_libraries(states harmonic)
target_link_libraries(states utils)
target_link_libraries(states common)

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(linear_algebra)

add_library(linear_algebra linear_algebra.F90 sparse_algebra.F90)
target_include_directories(linear_algebra PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against required modules.
target_link_libraries(linear_algebra precision)
target_link_libraries(linear_algebra macros)
target_link_libraries(linear_algebra io)
target_link_libraries(linear_algebra lapack)

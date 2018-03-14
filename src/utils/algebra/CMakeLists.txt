cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(algebra)

add_library(algebra mathematical_constants.f90 linear_algebra.f90 algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 eigenstuff.f90 integer_arrays.f90 group.f90 algebra.f90)

# Link against required modules.
add_subdirectory(${CMAKE_SOURCE_DIR}/utils/precision precision/)
target_link_libraries(algebra precision)

add_subdirectory(${CMAKE_SOURCE_DIR}/utils/logic logic/)
target_link_libraries(algebra logic)

add_subdirectory(${CMAKE_SOURCE_DIR}/utils/io io/)
target_link_libraries(algebra io)

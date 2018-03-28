cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(algebra)

add_library(algebra mathematical_constants.f90 linear_algebra.f90 algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 hermitian_eigenstuff.f90 orthonormal.f90 unitary_eigenstuff.f90 integer_arrays.f90 group.f90 algebra.f90)

# Link against required modules.
target_link_libraries(algebra precision)
target_link_libraries(algebra logic)
target_link_libraries(algebra io)

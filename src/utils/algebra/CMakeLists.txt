cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(algebra)

if (NOT DEFINED LINK_TO_LAPACK OR LINK_TO_LAPACK)
  add_library(algebra lapack_wrapper.f90 mathematical_constants.f90 linear_algebra.f90 algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 hermitian_eigenstuff.f90 orthonormal.f90 unitary_eigenstuff.f90 integer_arrays.f90 group.f90 integer_complex.f90 fraction_complex.f90 tests.f90 algebra.f90)

  # Link to LAPACK, if required.
  if (NOT DEFINED FIND_LAPACK OR FIND_LAPACK)
    find_package( LAPACK REQUIRED )
    target_link_libraries(algebra ${LAPACK_LIBRARIES})
  endif()
else()
  add_library(algebra lapack_dummy.f90 mathematical_constants.f90 linear_algebra.f90 algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 hermitian_eigenstuff.f90 orthonormal.f90 unitary_eigenstuff.f90 integer_arrays.f90 group.f90 integer_complex.f90 fraction_complex.f90 tests.f90 algebra.f90)
endif()

# Link against required modules.
target_link_libraries(algebra precision)
target_link_libraries(algebra abstract)
target_link_libraries(algebra io)

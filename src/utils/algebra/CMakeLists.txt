cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(algebra)

add_library(algebra algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 hermitian_eigenstuff.f90 orthonormal.f90 unitary_eigenstuff.f90 integer_arrays.f90 group.f90 integer_complex.f90 fraction_complex.f90 lanczos.f90 newton_raphson.f90 tests.f90 algebra.f90)

# Link against subsidiary libraries.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lapack lapack/)
target_link_libraries(algebra lapack)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/arpack arpack/)
target_link_libraries(algebra arpack)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mathematical_constants mathematical_constants/)
target_link_libraries(algebra mathematical_constants)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/linear_algebra linear_algebra/)
target_link_libraries(algebra linear_algebra)

# Link against required modules.
target_link_libraries(algebra precision)
target_link_libraries(algebra abstract)
target_link_libraries(algebra io)

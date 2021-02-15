cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(algebra)

add_library(algebra_modules algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 hermitian_eigenstuff.f90 orthonormal.f90 unitary_eigenstuff.f90 integer_arrays.f90 group.f90 integer_complex.f90 fraction_complex.f90 lanczos.f90 newton_raphson.f90 tests.f90 algebra.f90)
add_library(algebra_submodules ${DUMMY_FILE})
target_link_libraries(algebra_submodules algebra_modules)

# Link against subsidiary libraries.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lapack lapack/)
target_link_libraries(algebra_modules lapack_modules)
target_link_libraries(algebra_submodules lapack_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/arpack arpack/)
target_link_libraries(algebra_modules arpack_modules)
target_link_libraries(algebra_submodules arpack_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mathematical_constants mathematical_constants/)
target_link_libraries(algebra_modules mathematical_constants_modules)
target_link_libraries(algebra_submodules mathematical_constants_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/linear_algebra linear_algebra/)
target_link_libraries(algebra_modules linear_algebra_modules)
target_link_libraries(algebra_submodules linear_algebra_submodules)

# Link against required modules.
target_link_libraries(algebra_modules foundations_modules)
target_link_libraries(algebra_submodules foundations_submodules)
target_link_libraries(algebra_modules abstract_modules)
target_link_libraries(algebra_submodules abstract_submodules)
target_link_libraries(algebra_modules io_modules)
target_link_libraries(algebra_submodules io_submodules)

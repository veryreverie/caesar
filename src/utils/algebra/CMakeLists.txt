cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(algebra)

add_library(algebra_modules algebra_utils.f90 fraction.f90 phase.f90 fraction_algebra.f90 qr_decomposition.f90 hermitian_eigenstuff.f90 orthonormal.f90 unitary_eigenstuff.f90 integer_arrays.f90 group.f90 integer_complex.f90 fraction_complex.f90 lanczos.f90 newton_raphson.f90 tests.f90 algebra.f90)
add_library(algebra_submodules algebra_utils_submodule.f90 fraction_submodule.f90 phase_submodule.f90 fraction_algebra_submodule.f90 qr_decomposition_submodule.f90 hermitian_eigenstuff_submodule.f90 orthonormal_submodule.f90 unitary_eigenstuff_submodule.f90 integer_arrays_submodule.f90 group_submodule.f90 integer_complex_submodule.f90 fraction_complex_submodule.f90 lanczos_submodule.f90 newton_raphson_submodule.f90 tests_submodule.f90)
target_link_libraries(algebra_submodules algebra_modules)
add_dependencies(algebra_submodules algebra_modules)

# Link against subsidiary libraries.
foreach(dir IN ITEMS lapack arpack mathematical_constants linear_algebra)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(algebra_modules ${dir}_modules)
  target_link_libraries(algebra_submodules ${dir}_submodules)
  add_dependencies(algebra_modules ${dir}_modules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS foundations abstract io)
  target_link_libraries(algebra_modules ${dir}_modules)
  target_link_libraries(algebra_submodules ${dir}_submodules)
  add_dependencies(algebra_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs algebra_utils_test.pf fraction_test.pf phase_test.pf fraction_algebra_test.pf qr_decomposition_test.pf hermitian_eigenstuff_test.pf orthonormal_test.pf unitary_eigenstuff_test.pf integer_arrays_test.pf group_test.pf integer_complex_test.pf fraction_complex_test.pf lanczos_test.pf newton_raphson_test.pf tests_test.pf)
  add_pfunit_ctest(algebra_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(algebra_tests algebra_submodules)
endif()

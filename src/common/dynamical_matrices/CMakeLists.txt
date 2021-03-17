cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(dynamical_matrices)

# Normal mode modules.
add_library(dynamical_matrices_modules unique_directions.f90 min_images.f90 dynamical_matrix.f90 calculate_dynamical_matrices.f90 dynamical_matrices.f90)
add_library(dynamical_matrices_submodules unique_directions_submodule.f90 min_images_submodule.f90 dynamical_matrix_submodule.f90 calculate_dynamical_matrices_submodule.f90)
target_link_libraries(dynamical_matrices_submodules dynamical_matrices_modules)
add_dependencies(dynamical_matrices_submodules dynamical_matrices_modules)

# Link against required module.
foreach(dir IN ITEMS utils structure normal_mode)
  target_link_libraries(dynamical_matrices_modules ${dir}_modules)
  target_link_libraries(dynamical_matrices_submodules ${dir}_submodules)
  add_dependencies(dynamical_matrices_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs unique_directions_test.pf min_images_test.pf dynamical_matrix_test.pf calculate_dynamical_matrices_test.pf)
  add_pfunit_ctest(dynamical_matrices_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(dynamical_matrices_tests dynamical_matrices_submodules)
endif()

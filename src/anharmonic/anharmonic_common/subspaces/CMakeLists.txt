cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(subspaces)

add_library(subspaces_modules subspace_coupling.f90 subspace_combination.f90 subspaces.f90)
add_library(subspaces_submodules subspace_coupling_submodule.f90 subspace_combination_submodule.f90)
target_link_libraries(subspaces_submodules subspaces_modules)
add_dependencies(subspaces_submodules subspaces_modules)

# Link against required modules.
foreach(dir IN ITEMS common harmonic)
  target_link_libraries(subspaces_modules ${dir}_modules)
  target_link_libraries(subspaces_submodules ${dir}_submodules)
  add_dependencies(subspaces_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs subspace_coupling_test.pf subspace_combination_test.pf)
  add_pfunit_ctest(subspaces_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(subspaces_tests subspaces_submodules)
endif()

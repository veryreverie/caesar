cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(electronic_structure_common)

add_library(electronic_structure_common_modules kpoint_grid.f90 structure_file.f90 loto_splitting.f90 electronic_structure_common.f90)
add_library(electronic_structure_common_submodules kpoint_grid_submodule.f90 structure_file_submodule.f90 loto_splitting_submodule.f90)
target_link_libraries(electronic_structure_common_submodules electronic_structure_common_modules)

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode electronic_structure_data)
  target_link_libraries(electronic_structure_common_modules ${dir}_modules)
  target_link_libraries(electronic_structure_common_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs kpoint_grid_test.pf structure_file_test.pf loto_splitting_test.pf)
  add_pfunit_ctest(electronic_structure_common_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(electronic_structure_common_tests electronic_structure_common_submodules)
endif()

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(electronic_structure_data)

add_library(electronic_structure_data_modules linear_response.f90 electronic_structure_data.f90)
add_library(electronic_structure_data_submodules linear_response_submodule.f90 electronic_structure_data_submodule.f90)
target_link_libraries(electronic_structure_data_submodules electronic_structure_data_modules)

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode)
  target_link_libraries(electronic_structure_data_modules ${dir}_modules)
  target_link_libraries(electronic_structure_data_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs linear_response_test.pf electronic_structure_data_test.pf)
  add_pfunit_ctest(electronic_structure_data_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(electronic_structure_data_tests electronic_structure_data_submodules)
endif()

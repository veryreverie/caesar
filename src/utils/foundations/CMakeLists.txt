cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(foundations)

add_library(foundations_modules precision.f90 traits.f90 foundations.f90)
add_library(foundations_submodules ${DUMMY_FILE})
target_link_libraries(foundations_submodules foundations_modules)

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs precision_test.pf)
  add_pfunit_ctest(foundations_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(foundations_tests foundations_submodules)
endif()

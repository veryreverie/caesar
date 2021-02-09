cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(precision)

add_library(precision precision.f90)

# Enable testing.
include(${CMAKE_SOURCE_DIR}/utils/link_pfunit.cmake)
if(BUILD_TESTS)
  set(test_srcs precision_test.pf)
  add_pfunit_ctest(precision_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(precision_tests precision)
endif()

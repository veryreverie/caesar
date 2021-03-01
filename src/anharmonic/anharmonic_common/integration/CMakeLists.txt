cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(integration)

# Normal mode modules.
add_library(integration_modules integratable.f90 integrated.f90 integration.f90)
add_library(integration_submodules integrated_submodule.f90)
target_link_libraries(integration_submodules integration_modules)

# Link against required modules.
foreach(dir IN ITEMS common)
  target_link_libraries(integration_modules ${dir}_modules)
  target_link_libraries(integration_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs integrated_test.pf)
  add_pfunit_ctest(integration_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(integration_tests integration_submodules)
endif()

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(abstract)

# Compile the correct shared_counter file.
if (NOT DEFINED CORRECT_COUNTER_BUG OR NOT CORRECT_COUNTER_BUG)
  add_library(shared_counter shared_counter_no_bug.f90)
else()
  add_library(shared_counter shared_counter_bug.f90)
endif()

# Compile module files.
add_library(abstract_modules traits.f90 logic.f90 shared_counter.f90 abstract.f90)
target_link_libraries(abstract_modules shared_counter)

# Compile submodule files.
add_library(abstract logic_submodule.f90 shared_counter_submodule.f90)
target_link_libraries(abstract abstract_modules)

# Link against required modules.
target_link_libraries(abstract io_basic)

# Enable testing.
include(${CMAKE_SOURCE_DIR}/utils/link_pfunit.cmake)
if(BUILD_TESTS)
  set(test_srcs logic_test.pf shared_counter_test.pf)
  add_pfunit_ctest(abstract_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(abstract_tests abstract)
endif()

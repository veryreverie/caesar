cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(abstract)

# Compile the correct shared_counter file.
if (NOT DEFINED CORRECT_COUNTER_BUG OR NOT CORRECT_COUNTER_BUG)
  add_library(shared_counter shared_counter_no_bug.f90)
else()
  add_library(shared_counter shared_counter_bug.f90)
endif()

add_library(abstract_modules logic.f90 shared_counter.f90 abstract.f90)
add_library(abstract_submodules logic_submodule.f90 shared_counter_submodule.f90)
target_link_libraries(abstract_modules shared_counter)
target_link_libraries(abstract_submodules abstract_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations)
  target_link_libraries(abstract_modules ${dir}_modules)
  target_link_libraries(abstract_submodules ${dir}_submodules)
endforeach()

# Link submodules to required submodules later in the module build tree.
foreach(dir IN ITEMS io)
  target_link_libraries(abstract_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs logic_test.pf shared_counter_test.pf)
  add_pfunit_ctest(abstract_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(abstract_tests abstract_submodules)
endif()

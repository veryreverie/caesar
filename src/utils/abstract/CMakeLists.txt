cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(abstract)

# Compile the correct version of shared_counter_bugfix.f90
if (NOT DEFINED CORRECT_COUNTER_BUG OR NOT CORRECT_COUNTER_BUG)
  add_library(abstract traits.f90 logic.f90 shared_counter_no_bug.f90 shared_counter.f90 abstract.f90)
else()
  add_library(abstract traits.f90 logic.f90 shared_counter_bug.f90 shared_counter.f90 abstract.f90)
endif()

# Link against required modules.
target_link_libraries(abstract io_basic)

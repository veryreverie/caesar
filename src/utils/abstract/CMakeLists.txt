cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(abstract)

add_library(abstract traits.f90 logic.f90 shared_counter.f90 abstract.f90)

# Link against required modules.
target_link_libraries(abstract io_basic)

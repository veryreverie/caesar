cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(random)

add_library(random random.f90)

# Link against required modules.
target_link_libraries(random precision)
target_link_libraries(random abstract)
target_link_libraries(random io)

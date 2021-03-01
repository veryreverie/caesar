cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(linear_algebra)

add_library(linear_algebra_modules linear_algebra.F90 sparse_algebra.F90)
target_include_directories(linear_algebra_modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_library(linear_algebra_submodules ${DUMMY_FILE})
target_link_libraries(linear_algebra_submodules linear_algebra_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations macros io lapack)
  target_link_libraries(linear_algebra_modules ${dir}_modules)
  target_link_libraries(linear_algebra_submodules ${dir}_submodules)
endforeach()

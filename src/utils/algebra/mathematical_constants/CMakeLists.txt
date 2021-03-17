cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(mathematical_constants)

add_library(mathematical_constants_modules mathematical_constants.f90)
add_library(mathematical_constants_submodules ${DUMMY_FILE})
target_link_libraries(mathematical_constants_submodules mathematical_constants_modules)
add_dependencies(mathematical_constants_submodules mathematical_constants_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations)
  target_link_libraries(mathematical_constants_modules ${dir}_modules)
  target_link_libraries(mathematical_constants_submodules ${dir}_submodules)
  add_dependencies(mathematical_constants_modules ${dir}_modules)
endforeach()

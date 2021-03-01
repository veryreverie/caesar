cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(xyz)

add_library(xyz_modules xyz.f90)
add_library(xyz_submodules ${DUMMY_FILE})
target_link_libraries(xyz_submodules xyz_modules)

# Link against required modules.
foreach(dir IN ITEMS utils structure)
  target_link_libraries(xyz_modules ${dir}_modules)
  target_link_libraries(xyz_submodules ${dir}_submodules)
endforeach()

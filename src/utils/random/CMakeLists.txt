cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(random)

add_library(random_modules random.f90)
add_library(random_submodules random_submodule.f90)
target_link_libraries(random_submodules random_modules)
add_dependencies(random_submodules random_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations abstract io)
  target_link_libraries(random_modules ${dir}_modules)
  target_link_libraries(random_submodules ${dir}_submodules)
  add_dependencies(random_modules ${dir}_modules)
endforeach()

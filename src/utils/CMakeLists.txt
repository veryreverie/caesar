cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(utils)

# Utility modules.
add_library(utils_modules utils.f90)
add_library(utils_submodules ${DUMMY_FILE})
target_link_libraries(utils_submodules utils_modules)
add_dependencies(utils_submodules utils_modules)

# Link to subsidiary libraries.
foreach(dir IN ITEMS foundations macros abstract io algebra random arguments)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(utils_modules ${dir}_modules)
  target_link_libraries(utils_submodules ${dir}_submodules)
  add_dependencies(utils_modules ${dir}_modules)
endforeach()

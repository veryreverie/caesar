cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(quip)

add_library(quip_modules quip.f90)
add_library(quip_submodules ${DUMMY_FILE})
target_link_libraries(quip_submodules quip_modules)
add_dependencies(quip_submodules quip_modules)

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/quip_wrapper quip_wrapper/)
target_link_libraries(quip_modules quip_wrapper)
add_dependencies(quip_modules quip_wrapper)

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode electronic_structure_data)
  target_link_libraries(quip_modules ${dir}_modules)
  target_link_libraries(quip_submodules ${dir}_submodules)
  add_dependencies(quip_modules ${dir}_modules)
endforeach()

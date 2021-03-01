cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(common)

# Common modules.
add_library(common_modules common.f90)
add_library(common_submodules ${DUMMY_FILE})
target_link_libraries(common_submodules common_modules)

# Include subsidiary modules.
foreach(dir IN ITEMS structure normal_mode dynamical_matrices harmonic_stress electronic_structure observables)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(common_modules ${dir}_modules)
  target_link_libraries(common_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS utils)
  target_link_libraries(common_modules ${dir}_modules)
  target_link_libraries(common_submodules ${dir}_submodules)
endforeach()

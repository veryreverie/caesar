cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(castep)

add_library(castep_modules castep_wrapper.f90 phonon_file.f90 bands.f90 castep.f90)
add_library(castep_submodules ${DUMMY_FILE})
target_link_libraries(castep_submodules castep_modules)
add_dependencies(castep_submodules castep_modules)

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode dynamical_matrices electronic_structure_data electronic_structure_common)
  target_link_libraries(castep_modules ${dir}_modules)
  target_link_libraries(castep_submodules ${dir}_submodules)
  add_dependencies(castep_modules ${dir}_modules)
endforeach()

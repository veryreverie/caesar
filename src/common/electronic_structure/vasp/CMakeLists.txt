cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(vasp)

add_library(vasp_modules vasp_wrapper.f90 vasp.f90)
add_library(vasp_submodules ${DUMMY_FILE})
target_link_libraries(vasp_submodules vasp_modules)
add_dependencies(vasp_submodules vasp_modules)

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode dynamical_matrices electronic_structure_data electronic_structure_common)
  target_link_libraries(vasp_modules ${dir}_modules)
  target_link_libraries(vasp_submodules ${dir}_submodules)
  add_dependencies(vasp_modules ${dir}_modules)
endforeach()

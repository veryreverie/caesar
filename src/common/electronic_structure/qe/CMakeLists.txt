cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(qe)

add_library(qe_modules qe_fc.f90 qe_pw.f90 qe.f90)
add_library(qe_submodules ${DUMMY_FILE})
target_link_libraries(qe_submodules qe_modules)
add_dependencies(qe_submodules qe_modules)

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode dynamical_matrices electronic_structure_data electronic_structure_common)
  target_link_libraries(qe_modules ${dir}_modules)
  target_link_libraries(qe_submodules ${dir}_submodules)
  add_dependencies(qe_modules ${dir}_modules)
endforeach()

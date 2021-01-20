cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(atom)

add_library(atom basic_atom.f90 basic_structure.f90 atom_data.f90 atom.f90)

# Link against required modules.
target_link_libraries(atom utils)

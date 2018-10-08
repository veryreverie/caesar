cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(structure)

add_library(structure physical_constants.f90 basic_symmetry.f90 basic_structure.f90 qpoint.f90 generate_qpoints.f90 symmetry.f90 structure_data.f90 supercell.f90 structure.f90)

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/atom atom/)
target_link_libraries(structure atom)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/spglib spglib/)
target_link_libraries(structure spglib)

# Link against required modules.
target_link_libraries(structure utils)

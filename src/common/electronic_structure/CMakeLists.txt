cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(electronic_structure)

add_library(electronic_structure structure_file.f90 castep_wrapper.f90 qe_wrapper.f90 vasp_wrapper.f90 electronic_structure_file.f90 calculation_writer.f90 calculation_runner.f90 calculation_reader.f90 phonon_file.f90 electronic_structure.f90)

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/electronic_structure_data electronic_structure_data/)
target_link_libraries(electronic_structure electronic_structure_data)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/quip quip/)
target_link_libraries(electronic_structure quip)

# Link against required modules.
target_link_libraries(electronic_structure utils)
target_link_libraries(electronic_structure structure)
target_link_libraries(electronic_structure normal_mode)

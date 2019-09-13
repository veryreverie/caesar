cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(electronic_structure)

add_library(electronic_structure electronic_structure_file.f90 calculation_writer.f90 calculation_runner.f90 calculation_reader.f90 electronic_structure.f90)

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/electronic_structure_data electronic_structure_data/)
target_link_libraries(electronic_structure electronic_structure_data)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/electronic_structure_common electronic_structure_common/)
target_link_libraries(electronic_structure electronic_structure_common)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/quip quip/)
target_link_libraries(electronic_structure quip)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/castep castep/)
target_link_libraries(electronic_structure castep)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/qe qe/)
target_link_libraries(electronic_structure qe)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vasp vasp/)
target_link_libraries(electronic_structure vasp)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/xyz xyz/)
target_link_libraries(electronic_structure xyz)

# Link against required modules.
target_link_libraries(electronic_structure utils)
target_link_libraries(electronic_structure structure)
target_link_libraries(electronic_structure normal_mode)
target_link_libraries(electronic_structure dynamical_matrices)

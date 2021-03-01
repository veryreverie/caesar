cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(electronic_structure)

add_library(electronic_structure_modules electronic_structure_file.f90 calculation_writer.f90 calculation_runner.f90 calculation_reader.f90 electronic_structure.f90)
add_library(electronic_structure_submodules electronic_structure_file_submodule.f90 calculation_writer_submodule.f90 calculation_runner_submodule.f90 calculation_reader_submodule.f90)
target_link_libraries(electronic_structure_submodules electronic_structure_modules)

# Include subsidiary modules.
foreach(dir IN ITEMS electronic_structure_data electronic_structure_common quip castep qe vasp xyz)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(electronic_structure_modules ${dir}_modules)
  target_link_libraries(electronic_structure_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS utils structure normal_mode dynamical_matrices)
  target_link_libraries(electronic_structure_modules ${dir}_modules)
  target_link_libraries(electronic_structure_submodules ${dir}_submodules)
endforeach()

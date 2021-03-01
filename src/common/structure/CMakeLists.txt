cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(structure)

add_library(structure_modules physical_constants.f90 basic_symmetry.f90 qpoint.f90 generate_qpoints.f90 symmetry.f90 structure_data.f90 supercell.f90 structure.f90)
add_library(structure_submodules physical_constants_submodule.f90 basic_symmetry_submodule.f90 qpoint_submodule.f90 generate_qpoints_submodule.f90 symmetry_submodule.f90 structure_data_submodule.f90 supercell_submodule.f90)
target_link_libraries(structure_submodules structure_modules)

# Include subsidiary modules.
foreach(dir IN ITEMS atom spglib)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(structure_modules ${dir}_modules)
  target_link_libraries(structure_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS utils)
  target_link_libraries(structure_modules ${dir}_modules)
  target_link_libraries(structure_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs physical_constants_test.pf basic_symmetry_test.pf qpoint_test.pf generate_qpoints_test.pf symmetry_test.pf structure_data_test.pf supercell_test.pf)
  add_pfunit_ctest(structure_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(structure_tests structure_submodules)
endif()

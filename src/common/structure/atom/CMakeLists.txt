cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(atom)

add_library(atom_modules basic_atom.f90 basic_structure.f90 atom_data.f90 atom.f90)
add_library(atom_submodules basic_atom_submodule.f90 basic_structure_submodule.f90 atom_data_submodule.f90)
target_link_libraries(atom_submodules atom_modules)
add_dependencies(atom_submodules atom_modules)

# Link against required modules.
foreach(dir IN ITEMS utils)
  target_link_libraries(atom_modules ${dir}_modules)
  target_link_libraries(atom_submodules ${dir}_submodules)
  add_dependencies(atom_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs basic_atom_test.pf basic_structure_test.pf atom_data_test.pf)
  add_pfunit_ctest(atom_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(atom_tests atom_submodules)
endif()

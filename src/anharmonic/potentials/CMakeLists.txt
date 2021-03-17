cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(potentials)

add_library(potentials_modules potential_example.f90 potentials.f90)
add_library(potentials_submodules potential_example_submodule.f90 potential_data_types_submodule.f90 stress_data_types_submodule.f90 abstract_potential_types_submodule.f90)
target_link_libraries(potentials_submodules potentials_modules)
add_dependencies(potentials_submodules potentials_modules)
target_link_libraries(potentials_modules potentials_submodules)

# Include subsidiary modules.
foreach(dir IN ITEMS polynomial)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(potentials_modules ${dir}_modules)
  target_link_libraries(potentials_submodules ${dir}_submodules)
  add_dependencies(potentials_modules ${dir}_modules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS common anharmonic_common states)
  target_link_libraries(potentials_modules ${dir}_modules)
  target_link_libraries(potentials_submodules ${dir}_submodules)
  add_dependencies(potentials_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs potential_example_test.pf)
  add_pfunit_ctest(potentials_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(potentials_tests potentials_submodules)
endif()

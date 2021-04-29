cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(stars)

add_library(stars_modules qpoint_power.f90 qpoint_combination.f90 qpoint_star.f90 subspace_qpoint_stars.f90 stars.f90)
add_library(stars_submodules qpoint_power_submodule.f90 qpoint_combination_submodule.f90 qpoint_combinations_submodule.f90 qpoint_star_submodule.f90 qpoint_stars_submodule.f90 subspace_qpoint_stars_submodule.f90)
target_link_libraries(stars_submodules stars_modules)
add_dependencies(stars_submodules stars_modules)

# Link against required modules.
foreach(dir IN ITEMS common)
  target_link_libraries(stars_modules ${dir}_modules)
  target_link_libraries(stars_submodules ${dir}_submodules)
  add_dependencies(stars_modules ${dir}_modules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs qpoint_power_test.pf qpoint_combination_test.pf qpoint_combinations_test.pf qpoint_star_test.pf qpoint_stars_test.pf subspace_qpoint_stars_test.pf)
  add_pfunit_ctest(stars_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(stars_tests stars_submodules)
endif()

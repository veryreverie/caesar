cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(effective_harmonic)

add_library(effective_harmonic_modules harmonic_states.f90 harmonic_basis.f90 effective_harmonic.f90)
add_library(effective_harmonic_submodules harmonic_states_submodule.f90 harmonic_basis_submodule.f90)
target_link_libraries(effective_harmonic_submodules effective_harmonic_modules)
target_include_directories(effective_harmonic_modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link against other modules.
foreach(dir IN ITEMS common anharmonic_common)
  target_link_libraries(effective_harmonic_modules ${dir}_modules)
  target_link_libraries(effective_harmonic_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs harmonic_states_test.pf harmonic_basis_test.pf)
  add_pfunit_ctest(effective_harmonic_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(effective_harmonic_tests effective_harmonic_submodules)
endif()

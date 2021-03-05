cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(arguments)

add_library(arguments_modules keyword.f90 common_keywords.f90 help.f90 dictionary.f90 program_mode.f90 program_modes.f90 arguments.f90)
add_library(arguments_submodules keyword_submodule.f90 common_keywords_submodule.f90 help_submodule.f90 dictionary_submodule.f90 program_mode_submodule.f90 program_modes_submodule.f90)
target_link_libraries(arguments_submodules arguments_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations abstract io random)
  target_link_libraries(arguments_modules ${dir}_modules)
  target_link_libraries(arguments_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs keyword_test.pf common_keywords_test.pf help_test.pf dictionary_test.pf program_mode_test.pf program_modes_test.pf)
  add_pfunit_ctest(arguments_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(arguments_tests arguments_submodules)
endif()

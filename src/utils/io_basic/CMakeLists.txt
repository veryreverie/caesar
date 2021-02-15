cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(io_basic)

# Compile compiler-specific functionality.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific compiler_specific_gfortran.f90)
else()
  add_library(compiler_specific compiler_specific_default.f90)
endif()

# Compile module files.
add_library(io_basic_modules system.c terminal.f90 error.f90 string_base.f90 string.f90 print_settings.f90 print.f90 intrinsics.f90 io_utils.f90 token.f90 io_basic.f90)
target_link_libraries(io_basic_modules compiler_specific)

# Compile submodule files.
add_library(io_basic error_submodule.f90 string_base_submodule.f90 string_submodule.f90 print_settings_submodule.f90 print_submodule.f90 intrinsics_submodule.f90 io_utils_submodule.f90 token_submodule.f90)
target_link_libraries(io_basic io_basic_modules)

# Link against required modules.
target_link_libraries(io_basic precision)

# Enable testing.
include(${CMAKE_SOURCE_DIR}/utils/link_pfunit.cmake)
if(BUILD_TESTS)
  set(test_srcs string_base_test.pf string_test.pf print_settings_test.pf intrinsics_test.pf token_test.pf)
  add_pfunit_ctest(io_basic_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(io_basic_tests io_basic)
endif()

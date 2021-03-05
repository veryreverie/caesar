cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(io)

# Compile compiler-specific functionality.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific compiler_specific_gfortran.f90)
else()
  add_library(compiler_specific compiler_specific_default.f90)
endif()

add_definitions(-DEXECUTABLE_LOCATION="${CMAKE_BINARY_DIR}/bin" -DPYTHON_SCRIPTS_LOCATION="${CMAKE_BINARY_DIR}/python")

add_library(io_modules system.c terminal.f90 error.f90 string_base.f90 string.f90 print_settings.f90 print.f90 intrinsics.f90 c_string.f90 home_directory.f90 current_working_directory.f90 executable_locations.f90 io_utils.f90 token.f90 string_array.f90 string_writeable.f90 strings_writeable.f90 file.f90 ifile.f90 ofile_target.f90 ofile.f90 string_readable.f90 strings_readable.f90 stringable.f90 stringsable.f90 io.f90)
add_library(io_submodules error_submodule.f90 string_base_submodule.f90 string_submodule.f90 print_settings_submodule.f90 print_submodule.f90 intrinsics_submodule.f90 c_string_submodule.f90 home_directory_submodule.f90 current_working_directory_submodule.f90 executable_locations_submodule.F90 io_utils_submodule.f90 token_submodule.f90 string_array_submodule.f90 string_writeable_submodule.f90 strings_writeable_submodule.f90 file_submodule.f90 ifile_submodule.f90 ofile_target_submodule.f90 ofile_submodule.f90)
target_link_libraries(io_modules compiler_specific)
target_link_libraries(io_submodules io_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations abstract)
  target_link_libraries(io_modules ${dir}_modules)
  target_link_libraries(io_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs string_base_test.pf string_test.pf print_settings_test.pf intrinsics_test.pf token_test.pf string_array_test.pf string_writeable_test.pf strings_writeable_test.pf string_readable_test.pf strings_readable_test.pf stringable_test.pf stringsable_test.pf)
  add_pfunit_ctest(io_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(io_tests io_submodules)
endif()

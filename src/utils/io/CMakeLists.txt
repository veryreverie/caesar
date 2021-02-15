cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(io)

# Compile module files.
add_library(io_modules string_array.f90 string_writeable.f90 strings_writeable.f90 file.f90 ifile.f90 ofile_target.f90 ofile.f90 string_readable.f90 strings_readable.f90 stringable.f90 stringsable.f90 io.f90)

# Compile submodule files.
add_library(io string_array_submodule.f90 string_writeable_submodule.f90 strings_writeable_submodule.f90 file_submodule.f90 ifile_submodule.f90 ofile_target_submodule.f90 ofile_submodule.f90)
target_link_libraries(io io_modules)

# Link against required modules.
target_link_libraries(io_modules precision)
target_link_libraries(io_modules io_basic)
target_link_libraries(io_modules abstract)

# Enable testing.
include(${CMAKE_SOURCE_DIR}/utils/link_pfunit.cmake)
if(BUILD_TESTS)
  set(test_srcs string_array_test.pf string_writeable_test.pf strings_writeable_test.pf string_readable_test.pf strings_readable_test.pf stringable_test.pf stringsable_test.pf)
  add_pfunit_ctest(io_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(io_tests io)
endif()

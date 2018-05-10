cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(io)

add_library(io string_array.f90 string_writeable.f90 strings_writeable.f90 file.f90 ifile.f90 ofile_target.f90 ofile.f90 string_readable.f90 strings_readable.f90 stringable.f90 stringsable.f90 io.f90)

# Link against required modules.
target_link_libraries(io precision)
target_link_libraries(io io_basic)
target_link_libraries(io abstract)

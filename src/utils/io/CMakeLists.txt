cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(io)

add_library(io system.c terminal.f90 error.f90 string_base.f90 string.f90 print.f90 intrinsics.f90 stringable.f90 printable.f90 string_array.f90 io_main.f90 file.f90 ifile.f90 ofile.f90 io.f90)

# Compile compiler-specific functionality.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific compiler_specific_gfortran.f90)
else()
  add_library(compiler_specific compiler_specific_default.f90)
endif()
target_link_libraries(io compiler_specific)

# Link against required modules.
target_link_libraries(io precision)

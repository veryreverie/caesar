cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(utils)

# Compile compiler-specific functionality.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific compiler_specific_gfortran.f90)
else()
  add_library(compiler_specific compiler_specific_default.f90)
endif()

# Utility modules.
add_library(utils system.c spglib_wrapper.c constants.f90 terminal.f90 error.f90 string_base.f90 string.f90 stringable.f90 printable.f90 io.f90 file.f90 ifile.f90 ofile.f90 utils.f90 logic.f90 keyword.f90 help.f90 dictionary.f90 fraction.f90 fraction_algebra.f90 phase.f90 linear_algebra.f90 qr_decomposition.f90 eigenstuff.f90 integer_arrays.f90 group.f90 process_arguments.f90 caesar_modes.f90)
target_link_libraries(utils compiler_specific)

# Link to LAPACK.
find_package( LAPACK REQUIRED )
target_link_libraries(utils ${LAPACK_LIBRARIES})

# Link to spglib.
find_library(spglib_library_path symspg)
if("${spglib_libray_path}" STREQUAL "spglib_library_path-NOTFOUND")
  MESSAGE(FATAL_ERROR "libsymspg not found. Please ensure spglib is on LIB path.")
endif()
target_link_libraries(utils ${spglib_library_path})

find_path(spglib_include_path spglib.h)
if("${spglib_include_path}" STREQUAL "spglib_include_path-NOTFOUND")
  MESSAGE(FATAL_ERROR "spglib.h not found. Please ensure spglib is on PATH.")
endif()
target_include_directories(utils PUBLIC ${spglib_include_path})

target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

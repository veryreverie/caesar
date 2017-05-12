cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(utils)

# Compile compiler-specific functionality.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific compiler_specific_gfortran.f90)
else()
  add_library(compiler_specific compiler_specific_default.f90)
endif()

add_library(utils system.c spglib_wrapper.c constants.f90 string.f90 io.f90 utils.f90 help.f90 dictionary.f90 linear_algebra.f90 group.f90 qpoints.f90 moller_plesset.f90 structure.f90 dft_input_file.f90 dft_output_file.f90 bands.f90 normal_mode.f90)
target_link_libraries(utils compiler_specific)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


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
include_directories(${spglib_include_path})

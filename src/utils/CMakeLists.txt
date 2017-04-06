cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(utils)

add_library(utils system.c spglib_wrapper.c constants.f90 string.f90 io.f90 utils.f90 linear_algebra.f90 algebra.f90 group.f90 supercell.f90 kpoints.f90 moller_plesset.f90 structure.f90 dft_input_file.f90 dft_output_file.f90 structure_to_dft.f90 calculate_symmetry.f90 bands.f90 displacement_patterns.f90)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Compiler compiler-specific functionality.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific gfortran_specific.f90)
endif()
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
include_directories(${spglib_include_path})

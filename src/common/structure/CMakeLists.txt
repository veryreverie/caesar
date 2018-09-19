cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(structure)

if (NOT DEFINED LINK_TO_SPGLIB OR LINK_TO_SPGLIB)
  add_library(structure physical_constants.f90 basic_symmetry.f90 basic_structure.f90 atom.f90 qpoint.f90 generate_qpoints.f90 spglib_wrapper.c spglib_wrapper.f90 symmetry.f90 structure_main.f90 supercell.f90 structure.f90)

  # Link to spglib.
  find_library(spglib_library_path symspg)
  if("${spglib_libray_path}" STREQUAL "spglib_library_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "libsymspg not found. Please ensure spglib is on LIB path.")
  endif()
  target_link_libraries(structure ${spglib_library_path})

  find_path(spglib_include_path spglib.h)
  if("${spglib_include_path}" STREQUAL "spglib_include_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "spglib.h not found. Please ensure spglib is on PATH.")
  endif()
  target_include_directories(structure PUBLIC ${spglib_include_path})
else()
  add_library(structure physical_constants.f90 basic_symmetry.f90 basic_structure.f90 atom.f90 qpoint.f90 generate_qpoints.f90 spglib_dummy.f90 symmetry.f90 structure_main.f90 supercell.f90 structure.f90)
endif()

# Link against required modules.
target_link_libraries(structure utils)

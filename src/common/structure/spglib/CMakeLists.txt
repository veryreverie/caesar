cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(spglib)

if (NOT DEFINED LINK_TO_SPGLIB OR LINK_TO_SPGLIB)
  add_library(spglib spglib_symmetries.f90 spglib_wrapper.c spglib_wrapper.f90 spglib.f90)

  # Link to spglib.
  find_library(spglib_library_path symspg)
  if("${spglib_libray_path}" STREQUAL "spglib_library_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "libsymspg not found. Please ensure spglib is on LIB path. See README.txt for more details.")
  endif()
  if("${spglib_libray_path}" STREQUAL "NOTFOUND")
    MESSAGE(FATAL_ERROR "libsymspg not found. Please ensure spglib is on LIB path. See README.txt for more details.")
  endif()
  target_link_libraries(spglib ${spglib_library_path})

  find_path(spglib_include_path spglib.h)
  if("${spglib_include_path}" STREQUAL "spglib_include_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "spglib.h not found. Please ensure spglib is on PATH. See README.txt for more details.")
  endif()
  if("${spglib_include_path}" STREQUAL "NOTFOUND")
    MESSAGE(FATAL_ERROR "spglib.h not found. Please ensure spglib is on PATH. See README.txt for more details.")
  endif()
  target_include_directories(spglib PUBLIC ${spglib_include_path})
else()
  add_library(spglib spglib_symmetries.f90 spglib_dummy.f90 spglib.f90)
endif()

# Link against required modules.
target_link_libraries(spglib utils)
target_link_libraries(spglib atom)

# ======================================================================
# Provides various utilities for use with cmake.
# ======================================================================

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

# A macro for include guards in other files.
macro(cmake_include_guard name)
  get_property(already_included GLOBAL PROPERTY ${name}_included SET)
  if(${already_included})
    return()
  else()
  endif()
  set_property(GLOBAL PROPERTY ${name}_included true)
endmacro()

# Find and link pFunit if ENABLE_TESTS is true.
if(NOT DEFINED ENABLE_TESTS OR ENABLE_TESTS)
  enable_language(Fortran)
  enable_language(C)
  if(NOT "${PFUNIT_FOUND}")
    find_package(PFUNIT REQUIRED)
  endif()
  enable_testing()
  link_directories("${PFUNIT_INCLUDE_DIRS}/../lib")
  set(BUILD_TESTS true)
else()
  set(BUILD_TESTS false)
endif()

# Add dummy.f90 so that libraries containing no source files can be made.
if(NOT EXISTS "${CMAKE_BINARY_DIR}/dummy.f90")
  file(WRITE "${CMAKE_BINARY_DIR}/dummy.f90" "")
endif()
set(DUMMY_FILE "${CMAKE_BINARY_DIR}/dummy.f90")

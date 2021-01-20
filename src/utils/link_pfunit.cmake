# Including this file finds and links PFUNIT,
#    and enables CTest testing.
# This file must be included in every directory where tests are
#    wanted, and every parent directory of such directories,
#    recursively up to /src.
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

if (NOT DEFINED ENABLE_TESTS OR ENABLE_TESTS)
  find_package(PFUNIT REQUIRED)
  enable_testing()
  link_directories("${PFUNIT_INCLUDE_DIRS}/../lib")
  set(BUILD_TESTS true)
else()
  set(BUILD_TESTS false)
endif()

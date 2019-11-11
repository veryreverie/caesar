cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(arpack)

# Link to ARPACK, if required.
if (DEFINED LINK_TO_ARPACK AND LINK_TO_ARPACK)
  # Check that LAPACK has also been requested (ARPACK depends on BLAS).
  if (DEFINED LINK_TO_LAPACK AND NOT LINK_TO_LAPACK)
    MESSAGE(FATAL_ERROR "Linking to LAPACK has been disabled but linking to ARPACK has not. ARPACK requires LAPACK to run.")
  endif()
  add_library(arpack arpack_wrapper.f90)

  find_library(ARPACK_LIBRARY NAMES arpack PATHS "lib" "lib32" "lib64")
  if (ARPACK_LIBRARY)
    MESSAGE(STATUS "Found ARPACK library: ${ARPACK_LIBRARY}")
    target_link_libraries(arpack ${ARPACK_LIBRARY})
  else()
    MESSAGE(FATAL_ERROR "Could not find ARPACK library. Please ensure that libarpack.a is on LIB.")
  endif()
else()
  add_library(arpack arpack_dummy.f90)
endif()

# Link against required modules.
target_link_libraries(arpack precision)
target_link_libraries(arpack abstract)
target_link_libraries(arpack io)

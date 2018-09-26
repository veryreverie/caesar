cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(lapack)

# Link to LAPACK, if required.
if (NOT DEFINED LINK_TO_LAPACK OR LINK_TO_LAPACK)
  add_library(lapack lapack_wrapper.f90)

  # If LINK_TO_LAPACK is true but FIND_LAPACK is false,
  #    then LAPACK is requested but does not need to be linked,
  #    e.g. it is already bundled with the compiler.
  if (NOT DEFINED FIND_LAPACK OR FIND_LAPACK)
    find_package( LAPACK REQUIRED )
    target_link_libraries(lapack ${LAPACK_LIBRARIES})
  endif()
else()
  add_library(lapack lapack_dummy.f90)
endif()

# Link against required modules.
target_link_libraries(lapack precision)
target_link_libraries(lapack abstract)
target_link_libraries(lapack io)

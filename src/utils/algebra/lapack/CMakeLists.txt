cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(lapack)

# Link to LAPACK, if required.
if (NOT DEFINED LINK_TO_LAPACK OR LINK_TO_LAPACK)
  add_library(lapack_modules lapack_wrapper.f90)

  # If LINK_TO_LAPACK is true but FIND_LAPACK is false,
  #    then LAPACK is requested but does not need to be linked,
  #    e.g. it is already bundled with the compiler.
  if (NOT DEFINED FIND_LAPACK OR FIND_LAPACK)
    find_package( LAPACK REQUIRED )
    target_link_libraries(lapack_modules ${LAPACK_LIBRARIES})
  endif()
else()
  add_library(lapack_modules lapack_dummy.f90)
endif()
add_library(lapack_submodules ${DUMMY_FILE})
target_link_libraries(lapack_submodules lapack_modules)
add_dependencies(lapack_submodules lapack_modules)

# Link against required modules.
foreach(dir IN ITEMS foundations io)
  target_link_libraries(lapack_modules ${dir}_modules)
  target_link_libraries(lapack_submodules ${dir}_submodules)
  add_dependencies(lapack_modules ${dir}_modules)
endforeach()

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(utils)

# Utility modules.
add_library(utils utils.f90)

# Link to subsidiary libraries.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/precision precision/)
target_link_libraries(utils precision)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io io/)
target_link_libraries(utils io)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/logic logic/)
target_link_libraries(utils logic)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/arguments arguments/)
target_link_libraries(utils arguments)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/algebra algebra/)
target_link_libraries(utils algebra)

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

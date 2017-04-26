# CAESAR CmakeLists.txt
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(caesar)

enable_language(Fortran)
enable_language(C)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-g -O0 -W -Wall -Wextra -pedantic -fmax-errors=1 -fcheck=all -fbacktrace")
endif()

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_C_FLAGS "-g -O0 -W -Wall -Wextra -pedantic -fmax-errors=1 -std=c99")
endif()

# Add the main caesar executable as bin/caesar
project(caesar)
add_executable(caesar caesar.f90)

set_target_properties(caesar PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(caesar utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/harmonic harmonic/)
target_link_libraries(caesar harmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/quadratic quadratic/)
target_link_libraries(caesar quadratic)

add_subdirectory(${CMAKE_SOURCE_DIR}/testing ${CMAKE_BINARY_DIR}/testing/)
target_link_libraries(caesar testing)

# Compile the old code to bin/old/ for testing purposes.
add_subdirectory(${CMAKE_SOURCE_DIR}/old ${CMAKE_BINARY_DIR}/old/)

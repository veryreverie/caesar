# CAESAR CmakeLists.txt
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(caesar)

enable_language(Fortran)
enable_language(C)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/)

# N.B. to compile with a different compiler, call:
#    cmake /path/to/src -DCMAKE_Fortran_COMPILER:PATH=gfortran
#    cmake /path/to/src -DCMAKE_Fortran_COMPILER:PATH=ifort
#    cmake /path/to/src -DCMAKE_Fortran_COMPILER:PATH=nagfor

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-g -O0 -W -Wall -Wextra -pedantic -fmax-errors=1 -fcheck=all -fbacktrace -std=f2003 -fall-intrinsics -ffpe-trap=zero,overflow -fimplicit-none")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  set(CMAKE_Fortran_FLAGS "-f2003 -I ${CMAKE_BINARY_DIR}/mod/ -C=all -u")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "-stand f03")
else()
  message(FATAL_ERROR "Unknown compiler: ${CMAKE_Fortran_COMPILER_ID}")
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

add_subdirectory(${CMAKE_SOURCE_DIR}/castep castep/)
target_link_libraries(caesar castep)

add_subdirectory(${CMAKE_SOURCE_DIR}/harmonic harmonic/)
target_link_libraries(caesar harmonic)

#add_subdirectory(${CMAKE_SOURCE_DIR}/anharmonic anharmonic/)
#target_link_libraries(caesar anharmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/testing ${CMAKE_BINARY_DIR}/testing/)
target_link_libraries(caesar testing)

# Compile the old code to bin/old/ for testing purposes.
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_subdirectory(${CMAKE_SOURCE_DIR}/old ${CMAKE_BINARY_DIR}/old/)
endif()

# Copy the python scripts to bin/python/.
add_subdirectory(${CMAKE_SOURCE_DIR}/python ${CMAKE_BINARY_DIR}/python/)

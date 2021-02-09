# CAESAR CmakeLists.txt
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(caesar)

enable_language(Fortran)
enable_language(C)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/)

# N.B. to compile with a different compiler, call:
#    cmake /path/to/src -DCMAKE_Fortran_COMPILER:PATH=gfortran
#    cmake /path/to/src -DCMAKE_Fortran_COMPILER:PATH=ifort
#    cmake /path/to/src -DCMAKE_Fortran_COMPILER:PATH=nagfor
# N.B. this will require a compatible C compiler, set using the flag -DCMAKE_C_COMPILER in the same manner.

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  # Disabled warnings:
  #    maybe-uninitialized: erroneously triggers.
  #    unused-dummy-argument: triggers for inherited functions.
  set(CMAKE_Fortran_FLAGS "-g -O3 -march=native -W -Wall -Wextra -pedantic -fmax-errors=1  -std=f2008 -fall-intrinsics -ffpe-trap=zero,overflow -fimplicit-none -no-pie -Wno-maybe-uninitialized -Wno-unused-dummy-argument -fbacktrace -fno-omit-frame-pointer -ffree-line-length-none")
  if (CMAKE_Fortran_COMPILER_VERSION LESS 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -fcheck-array-temporaries")
  else()
    # Additional disabled warnings:
    #    function-elimination: this is desirable behaviour.
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all,no-array-temps -Wno-function-elimination")
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  set(CMAKE_Fortran_FLAGS "-f2008 -I ${CMAKE_BINARY_DIR}/mod/ -C=all -u -mtrace -nan -info -g -gline -colour")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  # Disabled warnings:
  #    5462: Global name too long: unnecessarily limits module names.
  #    7712: as unused-dummy-argument above.
  set(CMAKE_Fortran_FLAGS "-stand f08 -O2 -xHost -g -standard-realloc-lhs -standard-semantics -stand -traceback -warn all -diag-error-limit 1 -heap-arrays 0 -diag-disable=5462,7712 -implicitnone")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS "-C -c -g -traceback -O0")
else()
  message(FATAL_ERROR "Unknown compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_C_FLAGS "-g -O0 -W -Wall -Wextra -pedantic -fmax-errors=1 -std=c99")
endif()

# Add the main caesar executable as bin/caesar
project(caesar)
add_executable(caesar version.f90 caesar.f90)

set_target_properties(caesar PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Link against subsidiary modules.
add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(caesar utils)

add_subdirectory(${CMAKE_SOURCE_DIR}/common common/)
target_link_libraries(caesar common)

add_subdirectory(${CMAKE_SOURCE_DIR}/dft dft/)
target_link_libraries(caesar dft)

add_subdirectory(${CMAKE_SOURCE_DIR}/harmonic harmonic/)
target_link_libraries(caesar harmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/anharmonic anharmonic/)
target_link_libraries(caesar anharmonic)

add_subdirectory(${CMAKE_SOURCE_DIR}/testing ${CMAKE_BINARY_DIR}/testing/)
target_link_libraries(caesar testing)

# Compile the old code to bin/old/ for testing purposes.
#if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#  add_subdirectory(${CMAKE_SOURCE_DIR}/old ${CMAKE_BINARY_DIR}/old/)
#endif()

# Copy the python scripts to bin/python/.
add_subdirectory(${CMAKE_SOURCE_DIR}/python ${CMAKE_BINARY_DIR}/python/)

# Enable testing.
include(${CMAKE_SOURCE_DIR}/utils/link_pfunit.cmake)

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(utils)

add_library(utils system.c constants.f90 string.f90 file.f90 utils.f90 linear_algebra.f90 algebra.f90 group.f90 supercell.f90 kpoints.f90 moller_plesset.f90 structure.f90 dft_input_file.f90 dft_output_file.f90 structure_to_dft.f90 calculate_symmetry.f90 bands.f90 displacement_patterns.f90)
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  add_library(compiler_specific gfortran_specific.f90)
endif()
target_link_libraries(utils compiler_specific)

find_package( BLAS REQUIRED )
target_link_libraries(utils ${BLAS_LIBRARIES})

find_package( LAPACK REQUIRED )
target_link_libraries(utils ${LAPACK_LIBRARIES})


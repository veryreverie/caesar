cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(common)

# Link to QUIP, if PATH_TO_QUIP is defined.
if (DEFINED PATH_TO_QUIP)
  add_library(quip quip_wrapper.f90)

  target_link_libraries(quip utils)
  
  find_library(quip_library_path quip PATH ${PATH_TO_QUIP})
  if("${quip_library_path}" STREQUAL "quip_library_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "libquip not found. Please ensure that PATH_TO_QUIP is correct.")
  endif()
  target_link_libraries(quip ${quip_library_path})
  
  find_library(atoms_library_path atoms PATH ${PATH_TO_QUIP})
  if("${quip_library_path}" STREQUAL "atoms_library_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "libAtoms not found. Please ensure that PATH_TO_QUIP is correct.")
  endif()
  target_link_libraries(quip ${atoms_library_path})
  
  find_path(quip_include_path quip_unified_wrapper_module.mod PATH ${PATH_TO_QUIP})
  if("${quip_include_path}" STREQUAL "quip_include_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "quip_unified_wrapper_module.mod not found. Please ensure that PATH_TO_QUIP is correct.")
  endif()
  target_include_directories(quip PRIVATE ${quip_include_path})
  
  find_package( LAPACK REQUIRED )
  target_link_libraries(quip ${LAPACK_LIBRARIES})
else()
  add_library(quip quip_dummy.f90)
  target_link_libraries(quip utils)
endif()

# Common modules.
add_library(common qpoints.f90 atom.f90 basic_symmetry.f90 symmetry.f90 structure.f90 generate_qpoints.f90 input_file.f90 output_file.f90 bands.f90 normal_mode.f90 normal_mode_symmetry.f90)
target_link_libraries(common utils)
target_link_libraries(common quip)

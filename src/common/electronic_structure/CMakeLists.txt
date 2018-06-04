cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(electronic_structure)

if (DEFINED PATH_TO_QUIP)
  add_library(electronic_structure electronic_structure_data.f90 structure_file.f90 quip_wrapper.f90 castep_wrapper.f90 qe_wrapper.f90 vasp_wrapper.f90 electronic_structure_file.f90 mass_weighted_forces.f90 electronic_structure.f90)
else ()
  add_library(electronic_structure electronic_structure_data.f90 structure_file.f90 quip_dummy.f90 castep_wrapper.f90 qe_wrapper.f90 vasp_wrapper.f90 electronic_structure_file.f90 mass_weighted_forces.f90 electronic_structure.f90)
endif ()

# Link to QUIP, if PATH_TO_QUIP is defined.
if (DEFINED PATH_TO_QUIP)
  find_library(quip_library_path quip PATH ${PATH_TO_QUIP})
  if("${quip_library_path}" STREQUAL "quip_library_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "libquip not found. Please ensure that PATH_TO_QUIP is correct.")
  endif()
  target_link_libraries(electronic_structure ${quip_library_path})
  
  find_library(atoms_library_path atoms PATH ${PATH_TO_QUIP})
  if("${quip_library_path}" STREQUAL "atoms_library_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "libAtoms not found. Please ensure that PATH_TO_QUIP is correct.")
  endif()
  target_link_libraries(electronic_structure ${atoms_library_path})
  
  find_path(quip_include_path quip_unified_wrapper_module.mod PATH ${PATH_TO_QUIP})
  if("${quip_include_path}" STREQUAL "quip_include_path-NOTFOUND")
    MESSAGE(FATAL_ERROR "quip_unified_wrapper_module.mod not found. Please ensure that PATH_TO_QUIP is correct.")
  endif()
  target_include_directories(electronic_structure PRIVATE ${quip_include_path})
endif()

# Link against required modules.
target_link_libraries(electronic_structure utils)
target_link_libraries(electronic_structure structure)

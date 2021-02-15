cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(common)

# Common modules.
add_library(common bands.f90 common.f90)

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/structure structure/)
target_link_libraries(common structure)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/normal_mode normal_mode/)
target_link_libraries(common normal_mode)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/dynamical_matrices dynamical_matrices/)
target_link_libraries(common dynamical_matrices)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/harmonic_stress harmonic_stress/)
target_link_libraries(common harmonic_stress)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/electronic_structure electronic_structure/)
target_link_libraries(common electronic_structure)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/observables observables/)
target_link_libraries(common observables)

# Link against utils module.
target_link_libraries(common utils)

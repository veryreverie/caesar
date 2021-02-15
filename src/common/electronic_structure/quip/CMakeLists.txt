cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(quip)

add_library(quip quip.f90)

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/quip_wrapper quip_wrapper/)
target_link_libraries(quip quip_wrapper)

# Link against required modules.
target_link_libraries(quip quip_wrapper)
target_link_libraries(quip utils)
target_link_libraries(quip structure)
target_link_libraries(quip normal_mode)
target_link_libraries(quip electronic_structure_data)

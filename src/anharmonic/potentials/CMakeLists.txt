cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(potentials)

add_library(potentials potential_example.f90 potentials.f90)
target_include_directories(potentials PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/polynomial polynomial/)
target_link_libraries(potentials polynomial)

# Link against other modules.
target_link_libraries(potentials common)
target_link_libraries(potentials states)
target_link_libraries(potentials anharmonic_common)

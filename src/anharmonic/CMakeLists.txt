cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic)

add_library(anharmonic setup_anharmonic.f90 run_anharmonic.f90 calculate_potential.f90)
target_include_directories(anharmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/anharmonic_common anharmonic_common/)
target_link_libraries(anharmonic anharmonic_common)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/polynomial polynomial/)
target_link_libraries(anharmonic polynomial)

# Link against other modules.
target_link_libraries(anharmonic harmonic)
target_link_libraries(anharmonic utils)
target_link_libraries(anharmonic common)

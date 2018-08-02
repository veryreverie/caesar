cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic)

add_library(anharmonic setup_anharmonic.f90 run_anharmonic.f90 calculate_potential.f90 map_anharmonic_modes.f90 plot_effective_frequencies.f90 map_potential.f90 plot_potential_map.f90 calculate_states.f90)
target_include_directories(anharmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/states states/)
target_link_libraries(anharmonic states)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/anharmonic_common anharmonic_common/)
target_link_libraries(anharmonic anharmonic_common)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/polynomial polynomial/)
target_link_libraries(anharmonic polynomial)

# Link against other modules.
target_link_libraries(anharmonic harmonic)
target_link_libraries(anharmonic utils)
target_link_libraries(anharmonic common)

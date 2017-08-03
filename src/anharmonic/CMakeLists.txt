cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(anharmonic)

add_library(anharmonic sampling_point.f90 coupling.f90 setup_anharmonic.f90 run_anharmonic.f90 calculate_anharmonic.f90)
target_include_directories(anharmonic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(anharmonic utils)

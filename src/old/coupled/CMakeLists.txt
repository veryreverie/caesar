cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(coupled)

add_library(coupled spline2.f90 positions2.f90 coupled_atomic_units.f90 scf.f90 vscfCoupled.f90)
target_include_directories(coupled PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(${CMAKE_SOURCE_DIR}/utils utils/)
target_link_libraries(coupled utils)

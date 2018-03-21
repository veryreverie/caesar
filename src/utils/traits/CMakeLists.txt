cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(traits)

add_library(traits traits.f90)

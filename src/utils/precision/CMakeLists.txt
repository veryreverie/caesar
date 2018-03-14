cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(precision)

add_library(precision precision.f90)
target_include_directories(precision PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

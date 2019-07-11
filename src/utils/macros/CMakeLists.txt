# fortran_preprocessor Castep Driver CMakeLists.txt
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

add_library(macros macros.f90)
target_include_directories(macros PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

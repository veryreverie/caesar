# fortran_preprocessor Castep Driver CMakeLists.txt
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

add_library(macros_modules ${DUMMY_FILE})
add_library(macros_submodules ${DUMMY_FILE})
target_include_directories(macros_modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(macros_submodules macros_modules)
add_dependencies(macros_submodules macros_modules)

cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(arguments)

add_library(arguments keyword.f90 common_keywords.f90 help.f90 dictionary.f90 caesar_modes.f90 arguments.f90)

# Link against required modules.
target_link_libraries(arguments precision)
target_link_libraries(arguments abstract)
target_link_libraries(arguments io)
target_link_libraries(arguments random)

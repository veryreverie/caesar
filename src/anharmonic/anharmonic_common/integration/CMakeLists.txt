cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(integration)

# Normal mode modules.
add_library(integration integratable.f90 integrated.f90 integration.f90)

# Link against required module.
target_link_libraries(integration common)

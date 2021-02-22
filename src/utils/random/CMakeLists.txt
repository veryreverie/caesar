cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(random)

add_library(random_modules random.f90)
add_library(random_submodules random_submodule.f90)
target_link_libraries(random_submodules random_modules)

# Link against required modules.
target_link_libraries(random_modules foundations_modules)
target_link_libraries(random_submodules foundations_submodules)
target_link_libraries(random_modules abstract_modules)
target_link_libraries(random_submodules abstract_submodules)
target_link_libraries(random_modules io_modules)
target_link_libraries(random_submodules io_submodules)

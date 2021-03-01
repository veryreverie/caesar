cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(potentials)

add_library(potentials_modules potential_example.f90 potentials.f90)
add_library(potentials_submodules potential_example_submodule.f90 potentials_submodule.f90)
target_link_libraries(potentials_submodules potentials_modules)
target_include_directories(potentials_modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Include subsidiary modules.
foreach(dir IN ITEMS polynomial)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(potentials_modules ${dir}_modules)
  target_link_libraries(potentials_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS common anharmonic_common states)
  target_link_libraries(potentials_modules ${dir}_modules)
  target_link_libraries(potentials_submodules ${dir}_submodules)
endforeach()

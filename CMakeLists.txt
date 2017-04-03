# CAESAR CmakeLists.txt
cmake_minimum_required (version 2.8 fatal_error)
include(utils/cmake_include_guard)
include_guard(caesar)

enable_language(Fortran)
enable_language(C)
set (cmake_fortran_compiler gfortran)

project(caesar)
add_executable(caesar caesar.f90)

add_subdirectory(utils utils/)
target_link_libraries(caesar utils)

add_subdirectory(harmonic harmonic/)
target_link_libraries(caesar harmonic)

add_subdirectory(quadratic quadratic/)
target_link_libraries(caesar quadratic)

add_subdirectory(testing testing/)
target_link_libraries(caesar testing)

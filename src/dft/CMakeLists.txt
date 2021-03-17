cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(dft)

add_library(dft_modules castep_output_file.f90 converge_cutoff_and_kpoints.f90 plot_cutoff_and_kpoints.f90 converge_harmonic_frequencies.f90 plot_harmonic_convergence.f90 dft.f90)
add_library(dft_submodules castep_output_file_submodule.f90 converge_cutoff_and_kpoints_submodule.f90 plot_cutoff_and_kpoints_submodule.f90 converge_harmonic_frequencies_submodule.f90 plot_harmonic_convergence_submodule.f90 dft_submodule.f90)
target_link_libraries(dft_submodules dft_modules)
add_dependencies(dft_submodules dft_modules)

# Link against required modules.
foreach(dir IN ITEMS common harmonic)
  target_link_libraries(dft_modules ${dir}_modules)
  target_link_libraries(dft_submodules ${dir}_submodules)
  add_dependencies(dft_modules ${dir}_modules)
endforeach()

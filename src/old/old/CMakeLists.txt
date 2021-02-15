# Compiles the old caesar code, for testing purposes.
# Nothing in this directory links into the main caesar executable.
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(old)

find_package( LAPACK REQUIRED )

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-g -O0")
endif()

# Add harmonic executables.
foreach(exe combine_forces compare_kpoints construct_finite_displacement construct_matrix_force_cnsts construct_supercell convert_forces_from_Rybohr_to_eVang equilibrium_frac generate_kgrid generate_supercell_kpoint_mesh_qe generate_supercells)
  add_executable(${exe} harmonic/programs/${exe}.f90)
endforeach(exe)

foreach(exe fourier_interpolation lte lte_lower)
  add_executable(${exe} harmonic/programs/${exe}.f90)
  target_link_libraries(${exe} ${LAPACK_LIBRARIES})
endforeach(exe)

# Copy harmonic shell scripts.
foreach(shell clear_all convert_harmonic dyn_mats fetch_forces_castep fetch_forces_qe hartree_to_eV lte_harmonic rutgers_run_harmonic setup_harmonic structure_to_castep structure_to_qe structure_to_vasp tcm_cluster_run_harmonic)
  file(COPY harmonic/scripts/${shell}.sh DESTINATION "${CMAKE_BINARY_DIR}/old")
endforeach(shell)

# Add quadratic executables.
foreach(exe band_folding calculate_anharmonic calculate_bs calculate_gap generate_amplitudes generate_quadratic_configurations generate_sc_path quadratic_spline)
  add_executable(${exe} quadratic/programs/${exe}.f90)
endforeach(exe)

foreach(exe vscf_1d)
  add_executable(${exe} quadratic/programs/${exe}.f90)
  target_link_libraries(${exe} ${LAPACK_LIBRARIES})
endforeach(exe)

# Copy quadratic shell scripts.
foreach(shell anharmonic bs_quadratic convert_quadratic eigenval_castep_to_bands eigenval_vasp_to_bands setup_quadratic tcm_cleanup_anharmonic tcm_cleanup_bs tcm_cluster_run_quadratic)
  file(COPY quadratic/scripts/${shell}.sh DESTINATION "${CMAKE_BINARY_DIR}/old")
endforeach(shell)

# Copies python scripts to the binary directory.
# Nothing in this directory links into the main caesar executable.
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/cmake_include_guard.cmake)
cmake_include_guard(python)

# Copy python scripts.
foreach(script utils.py plot_cutoff_and_kpoints.py plot_normal_modes.py plot_dos_and_dispersion.py plot_thermodynamic_variables.py plot_harmonic_qpoint_convergence.py plot_modes.py plot_potential_map.py plot_vscf_modes.py plot_vscf_states.py plot_vscf_convergence.py plot_harmonic_convergence.py)
  set (input ${CMAKE_CURRENT_SOURCE_DIR}/${script})
  set (output ${CMAKE_CURRENT_BINARY_DIR}/${script})
  add_custom_command(
    COMMENT "Updating Python script '${script}'."
    DEPENDS ${input}
    OUTPUT ${output}
    COMMAND ${CMAKE_COMMAND} -E copy ${input} ${output}
  )
  add_custom_target(${script} DEPENDS ${output})
  add_dependencies(caesar ${script})
endforeach(script)

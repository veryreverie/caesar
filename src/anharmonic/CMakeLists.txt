cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(anharmonic)

add_library(anharmonic_modules generate_subspace_potentials.f90 mode_map.f90 vscf.f90 interpolation.f90 setup_anharmonic.f90 run_anharmonic.f90 calculate_potential.f90 map_modes.f90 plot_modes.f90 map_potential.f90 plot_potential_map.f90 map_vscf_modes.f90 plot_vscf_modes.f90 calculate_anharmonic_observables.f90 plot_vscf_states.f90 plot_vscf_convergence.f90 anharmonic.f90)
add_library(anharmonic_submodules generate_subspace_potentials_submodule.f90 mode_map_submodule.f90 vscf_submodule.f90 interpolation_submodule.f90 setup_anharmonic_submodule.f90 run_anharmonic_submodule.f90 calculate_potential_submodule.f90 map_modes_submodule.f90 plot_modes_submodule.f90 map_potential_submodule.f90 plot_potential_map_submodule.f90 map_vscf_modes_submodule.f90 plot_vscf_modes_submodule.f90 calculate_anharmonic_observables_submodule.f90 plot_vscf_states_submodule.f90 plot_vscf_convergence_submodule.f90 anharmonic_submodule.f90)
target_link_libraries(anharmonic_submodules anharmonic_modules)

# Include subsidiary modules.
foreach(dir IN ITEMS anharmonic_common states potentials)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir}/)
  target_link_libraries(anharmonic_modules ${dir}_modules)
  target_link_libraries(anharmonic_submodules ${dir}_submodules)
endforeach()

# Link against required modules.
foreach(dir IN ITEMS common harmonic)
  target_link_libraries(anharmonic_modules ${dir}_modules)
  target_link_libraries(anharmonic_submodules ${dir}_submodules)
endforeach()

# Enable testing.
if(BUILD_TESTS)
  set(test_srcs generate_subspace_potentials_test.pf mode_map_test.pf vscf_test.pf interpolation_test.pf setup_anharmonic_test.pf run_anharmonic_test.pf calculate_potential_test.pf map_modes_test.pf plot_modes_test.pf map_potential_test.pf plot_potential_map_test.pf map_vscf_modes_test.pf plot_vscf_modes_test.pf calculate_anharmonic_observables_test.pf plot_vscf_states_test.pf plot_vscf_convergence_test.pf)
  add_pfunit_ctest(anharmonic_tests TEST_SOURCES ${test_srcs} LINK_LIBRARIES funit ${PFUNIT_LIBRARIES})
  target_link_libraries(anharmonic_tests anharmonic_submodules)
endif()

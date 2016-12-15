program caesar
  use utils,   only : command_line_args, i2s
  use process, only : ProcessResult, system_process
  use file_io, only : open_write_file, open_read_file
  
  ! use harmonic modules
  use combine_forces_module,        only : combine_forces
  use compare_kpoints_module,       only : compare_kpoints
  use construct_finite_displacement_module,&
    &only : construct_finite_displacement
  use construct_matrix_force_cnsts_module,&
    &only : construct_matrix_force_cnsts
  use construct_supercell_module,   only : construct_supercell
  use convert_forces_from_Rybohr_to_eVang_module,&
    &only : convert_forces_from_Rybohr_to_eVang
  use equilibrium_frac_module,      only : equilibrium_frac
  use fourier_interpolation_module, only : fourier_interpolation
  use generate_kgrid_module,        only : generate_kgrid
  use generate_supercell_kpoint_mesh_qe_module,&
    &only : generate_supercell_kpoint_mesh_qe
  use generate_supercells_module,   only : generate_supercells
  use lte_module,                   only : lte
  use lte_lower_module,             only : lte_lower
  
  ! use quadratic modules
  use anharmonic_module,           only : anharmonic
  use band_folding_module,         only : band_folding
  use calculate_anharmonic_module, only : calculate_anharmonic
  use calculate_bs_module,         only : calculate_bs
  use calculate_gap_module,        only : calculate_gap
  use generate_amplitudes_module,  only : generate_amplitudes
  use generate_quadratic_configurations_module,&
    &only : generate_quadratic_configurations
  use generate_sc_path_module,     only : generate_sc_path
  use quadratic_spline_module,     only : quadratic_spline
  use vscf_1d_module,              only : vscf_1d
  
  implicit none
  
  integer                    :: i             ! loop index
  character(32), allocatable :: args(:)       ! command line arguments
  character(:),  allocatable :: arg           ! command line argument
  integer                    :: return_status ! system() status

  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 0) then
    write(*,*) 'No arguments given. For help, call caesar -h'
    deallocate(args)
    stop
  else
    allocate(character(len=len(trim(args(1)))) :: arg)
    arg = trim(args(1))
  endif
  
  if (arg == '-h' .or. arg == '--help') then
    write(*,*) 'caesar [-h] [option]'
    write(*,*) ''
    write(*,*) '-h :'
    write(*,*) '  Displays this help text'
    write(*,*) ''
    write(*,*) 'option : utilities :'
    write(*,*) '  hartree_to_eV :'
    write(*,*) '    Provides a Hartree to eV calculator'
    write(*,*) ''
    write(*,*) 'option : harmonic calculations :'
    write(*,*) '  setup_harmonic :'
    write(*,*) '    Sets up calculation'
    write(*,*) '  convert_harmonic :'
    write(*,*) '    Converts calculation to specific code'
    write(*,*) '    Choices are castep, vasp and quantum espresso'
    write(*,*) '    Should be called after setup_harmonic'
    write(*,*) '  tcm_cluster_run_harmonic :'
    write(*,*) '    Runs calculation on the TCM cluster'
    write(*,*) '    Should be called after convert_harmonic'
    write(*,*) '  rutgers_run_harmonic :'
    write(*,*) '    Runs calculations'
    write(*,*) '    Should be called after convert_harmonic'
    write(*,*) '  lte_harmonic :'
    write(*,*) '    Runs harmonic calculations'
    write(*,*) '    Should be run after one of the run_harmonic options'
    write(*,*) '  clear_all :'
    write(*,*) '    Deletes all temporary files and folders'
    write(*,*) ''
    write(*,*) 'option : quadratic calculations :'
    write(*,*) '  [quadratic help text yet to be written]'
  ! test
  elseif (arg == 'test') then
    i=1
    do return_status=1,10
      write(*,*) 'Testing: no.loops = '//trim(i2s(return_status))
      do i=i,i+return_status-1
        write(*,*) 'i = '//trim(i2s(i))
      enddo
    enddo
  ! Wrappers for Fortran 
  elseif (arg == 'band_folding') then
    call band_folding()
  elseif (arg == 'calculate_anharmonic') then
    call calculate_anharmonic()
  elseif (arg == 'calculate_bs') then
    call calculate_bs()
  elseif (arg == 'calculate_gap') then
    call calculate_gap()
  elseif (arg == 'combine_forces') then
    call combine_forces()
  elseif (arg == 'compare_kpoints') then
    call compare_kpoints()
  elseif (arg == 'construct_finite_displacement') then
    call construct_finite_displacement(args(2:))
  elseif (arg == 'construct_matrix_force_cnsts') then
    call construct_matrix_force_cnsts(args(2:))
  elseif (arg == 'construct_supercell') then
    call construct_supercell()
  elseif (arg == 'convert_forces_from_Rybohr_to_eVang') then
    call convert_forces_from_Rybohr_to_eVang()
  elseif (arg == 'equilibrium_frac') then
    call equilibrium_frac()
  elseif (arg == 'fourier_interpolation') then
    call fourier_interpolation()
  elseif (arg == 'generate_amplitudes') then
    call generate_amplitudes()
  elseif (arg == 'generate_kgrid') then
    call generate_kgrid(args(2:))
  elseif (arg == 'generate_quadratic_configurations') then
    call generate_quadratic_configurations(args(2:))
  elseif (arg == 'generate_sc_path') then
    call generate_sc_path()
  elseif (arg == 'generate_supercell_kpoint_mesh_qe') then
    call generate_supercell_kpoint_mesh_qe()
  elseif (arg == 'generate_supercells') then
    call generate_supercells()
  elseif (arg == 'lte') then
    call lte()
  elseif (arg == 'lte_lower') then
    call lte_lower()
  elseif (arg == 'quadratic_spline') then
    call quadratic_spline()
  elseif (arg == 'vscf_1d') then
    call vscf_1d()
  ! wrappers for shell scripts
  elseif (arg == 'anharmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'bs_quadratic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'clear_all') then
    return_status = system(arg//'.sh')
  elseif (arg == 'convert_harmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'convert_quadratic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'dyn_mats') then
    return_status = system(arg//'.sh')
  elseif (arg == 'eigenval_castep_to_bands') then
    return_status = system(arg//'.sh')
  elseif (arg == 'eigenval_vasp_to_bands') then
    return_status = system(arg//'.sh')
  elseif (arg == 'fetch_forces_castep') then
    return_status = system(arg//'.sh')
  elseif (arg == 'fetch_forces_qe') then
    return_status = system(arg//'.sh')
  elseif (arg == 'hartree_to_eV') then
    return_status = system(arg//'.sh')
  elseif (arg == 'lte_harmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'rutgers_run_harmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'setup_harmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'setup_quadratic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'structure_to_castep') then
    return_status = system(arg//'.sh')
  elseif (arg == 'structure_to_qe') then
    return_status = system(arg//'.sh')
  elseif (arg == 'structure_to_vasp') then
    return_status = system(arg//'.sh')
  elseif (arg == 'tcm_cleanup_anharmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'tcm_cleanup_bs') then
    return_status = system(arg//'.sh')
  elseif (arg == 'tcm_cluster_run_harmonic') then
    return_status = system(arg//'.sh')
  elseif (arg == 'tcm_cluster_run_quadratic') then
    return_status = system(arg//'.sh')
  ! wrappers for python scripts
  elseif (arg == 'get_kpoints') then
    return_status = system(arg//'.py')
  ! unrecognised argument
  else
    write(*,*) 'Unrecognised argument : '//arg
  endif
  
  deallocate(arg)
  deallocate(args)

end program

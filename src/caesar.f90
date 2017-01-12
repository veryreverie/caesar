program caesar
  ! use utility modules
  use constants, only : dp
  use utils,     only : command_line_args, i2s
  use file_io,   only : open_write_file, open_read_file
  
  ! use class modules
  use string_module
  
  ! use common modules
  use rundft_module,        only : rundft
  
  ! use harmonic modules
  use combine_forces_module,        only : combine_forces
  use compare_kpoints_module,       only : compare_kpoints
  use construct_finite_displacement_module,&
    &only : construct_finite_displacement
  use construct_matrix_force_cnsts_module,&
    &only : construct_matrix_force_cnsts
  use construct_supercell_module,   only : construct_supercell
  use equilibrium_frac_module,      only : equilibrium_frac
  use fourier_interpolation_module, only : fourier_interpolation
  use generate_kgrid_module,        only : generate_kgrid
  use generate_supercell_kpoint_mesh_qe_module,&
    &only : generate_supercell_kpoint_mesh_qe
  use generate_supercells_module,   only : generate_supercells
  use lte_module,                   only : lte
  use hartree_to_eV_module, only : hartree_to_eV
  
  ! use quadratic modules
  use anharmonic_module,           only : anharmonic
  use band_folding_module,         only : band_folding
  use calculate_bs_module,         only : calculate_bs
  use calculate_gap_module,        only : calculate_gap
  use generate_quadratic_configurations_module,&
    &only : generate_quadratic_configurations
  use generate_sc_path_module,     only : generate_sc_path
  
  implicit none
  
  integer                      :: i             ! loop index
  character(100),  allocatable :: args(:)       ! command line arguments
  character(:),    allocatable :: arg           ! command line argument
  character(1000), allocatable :: argstring     ! command line arguments
  integer                      :: return_status ! system() status
  
  ! lte variables
  real(dp) :: tol,tol2,delta
  
  ! testing variables
  type(String) :: temp_string
  type(String) :: temp_string2
  type(String) :: temp_string3
  
  ! read in command line arguments
  args = command_line_args()
  
  argstring = ''
  do i=2,size(args)
    argstring = trim(argstring)//' '//args(i)
  enddo
  
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
    write(*,*) '  lte_harmonic :'
    write(*,*) '    Runs harmonic calculations'
    write(*,*) '    Should be run after one of the run_harmonic options'
    write(*,*) '  clear_all :'
    write(*,*) '    Deletes all temporary files and folders'
    write(*,*) ''
    write(*,*) 'option : quadratic calculations :'
    write(*,*) '  setup_quadratic :'
    write(*,*) '    Sets up quadratic calculation'
    write(*,*) '    Should be run after lte_harmonic'
    write(*,*) '  convert_quadratic :'
    write(*,*) '    Converts calculation to specific code'
    write(*,*) '    Choices are castep, vasp and quantum espresso'
    write(*,*) '    Should be called after setup_quadratic'
    write(*,*) '  tcm_cluster_run_quadratic :'
    write(*,*) '    Runs calculation on the TCM cluster'
    write(*,*) '    Should be called after convert_quadratic'
    write(*,*) '  tcm_cleanup_anharmonic :'
    write(*,*) '    Collates energies from anharmonic calculations'
    write(*,*) '    Should be called after tcm_cluster_run_quadratic'
    write(*,*) '  anharmonic :'
    write(*,*) '    Runs anharmonic calculations'
    write(*,*) '    Should be called after tcm_cleanup_anharmonic'
    write(*,*) '  tcm_cleanup_bs'
    write(*,*) '    Collates bands from anharmonic calculations'
    write(*,*) '    Should be called after tcm_cluster_run_quadratic'
    write(*,*) '  bs_quadratic :'
    write(*,*) '    Runs band structure calculations'
    write(*,*) '    Should be called after tcm_cleanup_bs'
    write(*,*) '  get_kpoints :'
    write(*,*) '    [Help text pending]'
    write(*,*) '  eigenval_vasp_to_bands'
    write(*,*) '    [Help text pending]'
    write(*,*) '  calculate_gap'
    write(*,*) '    [Help text pending]'
  ! test
  elseif (arg == 'test') then
    temp_string = 'test'
    write(*,*) char(temp_string)
    temp_string = 'test2'
    write(*,*) char(temp_string)
    temp_string2 = temp_string
    write(*,*) char(temp_string2)
    temp_string2 = 12
    write(*,*) char(temp_string2)
    temp_string3 = temp_string//' '//temp_string2//' '//143
    write(*,*) char(temp_string3)
    call drop(temp_string)
  ! Wrappers for Fortran 
  elseif (arg == 'band_folding') then
    call band_folding(args(2:))
  elseif (arg == 'calculate_bs') then
    call calculate_bs(args(2:))
  elseif (arg == 'calculate_gap') then
    call calculate_gap()
  elseif (arg == 'combine_forces') then
    call combine_forces(args(2:))
  elseif (arg == 'compare_kpoints') then
    call compare_kpoints(args(2:))
  elseif (arg == 'construct_finite_displacement') then
    call construct_finite_displacement(args(2:))
  elseif (arg == 'construct_matrix_force_cnsts') then
    call construct_matrix_force_cnsts(args(2:))
  elseif (arg == 'construct_supercell') then
    call construct_supercell(args(2:))
  elseif (arg == 'equilibrium_frac') then
    call equilibrium_frac(args(2:))
  elseif (arg == 'fourier_interpolation') then
    call fourier_interpolation(args(2),args(3),args(4),args(5),args(6),args(7),&
      & args(8),args(9),args(10),args(11),args(12),args(13))
  elseif (arg == 'generate_kgrid') then
    call generate_kgrid(args(2:))
  elseif (arg == 'generate_quadratic_configurations') then
    call generate_quadratic_configurations(args(2:))
  elseif (arg == 'generate_sc_path') then
    call generate_sc_path(args(2:))
  elseif (arg == 'generate_supercell_kpoint_mesh_qe') then
    call generate_supercell_kpoint_mesh_qe(args(2:))
  elseif (arg == 'generate_supercells') then
    call generate_supercells(args(2:))
  elseif (arg == 'hartree_to_eV') then
    call hartree_to_eV()
  elseif (arg == 'lte') then
    read(args(2),*) tol
    read(args(3),*) tol2
    read(args(4),*) delta
    call lte(tol,tol2,delta,args(5),args(6),args(7),args(8),args(9),args(10),&
      & args(11),args(12),args(13),args(14),args(15),args(16),args(17),      &
      & args(18),args(19))
  elseif (arg == 'rundft') then
    call rundft(args(2:))
  ! wrappers for shell scripts
  elseif (arg == 'anharmonic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'bs_quadratic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'clear_all') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'convert_harmonic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'convert_quadratic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'eigenval_castep_to_bands') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'eigenval_vasp_to_bands') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'fetch_forces_castep') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'fetch_forces_qe') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'lte_harmonic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'setup_harmonic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'setup_quadratic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'structure_to_castep') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'structure_to_qe') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'structure_to_vasp') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'tcm_cleanup_anharmonic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'tcm_cleanup_bs') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'tcm_cluster_run_harmonic') then
    return_status = system(arg//'.sh '//trim(argstring))
  elseif (arg == 'tcm_cluster_run_quadratic') then
    return_status = system(arg//'.sh '//trim(argstring))
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

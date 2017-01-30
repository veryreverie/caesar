program caesar
  ! use utility modules
  use constants, only : dp
  use utils,     only : command_line_args, i2s
  use file_module,   only : open_write_file, open_read_file
  
  ! use class modules
  use string_module
  
  ! use common modules
  use rundft_module
  use structure_to_dft_module
  use calculate_symmetry_helper_module
  
  ! use harmonic modules
  use combine_forces_module
  use compare_kpoints_module
  use construct_finite_displacement_module
  use construct_matrix_force_cnsts_module
  use construct_supercell_module
  use equilibrium_frac_module
  use fourier_interpolation_module
  use generate_kgrid_module
  use generate_supercells_module
  use lte_module
  use fetch_forces_module
  
  ! use quadratic modules
  use setup_quadratic_module
  use anharmonic_module
  use bs_quadratic_module
  
  ! use misc modules
  use calculate_gap_module
  use hartree_to_eV_module
  
  implicit none
  
  integer                   :: i             ! loop index
  type(String), allocatable :: args(:)       ! command line arguments
  type(String)              :: argstring     ! command line arguments
  
  ! lte variables
  real(dp) :: tol,tol2,delta
  
  ! testing variables
  type(String) :: temp_string
  type(String) :: temp_string2
  type(String) :: temp_string3
  real(dp)     :: temp_real
  real(dp)     :: temp_real2
  
  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 0) then
    write(*,*) 'No arguments given. For help, call caesar -h'
    stop
  endif
  
  argstring = ''
  do i=2,size(args)
    argstring = argstring//' '//args(i)
  enddo
  
  if (args(1) == '-h' .or. args(1) == '--help') then
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
    write(*,*) '    Converts calculation to specific code'
    write(*,*) '    Choices are castep, vasp and quantum espresso'
    write(*,*) '  tcm_cluster_run_harmonic :'
    write(*,*) '    Runs calculation on the TCM cluster'
    write(*,*) '    Should be called after setup_harmonic'
    write(*,*) '  lte_harmonic :'
    write(*,*) '    Runs harmonic calculations'
    write(*,*) '    Should be run after tcm_cluster_run_harmonic'
    write(*,*) '  clear_all :'
    write(*,*) '    Deletes all temporary files and folders'
    write(*,*) ''
    write(*,*) 'option : quadratic calculations :'
    write(*,*) '  setup_quadratic :'
    write(*,*) '    Sets up quadratic calculation'
    write(*,*) '    Should be run after lte_harmonic'
    write(*,*) '    Converts calculation to specific code'
    write(*,*) '    Choices are castep, vasp and quantum espresso'
    write(*,*) '    Should be called after setup_quadratic'
    write(*,*) '  tcm_cluster_run_quadratic :'
    write(*,*) '    Runs calculation on the TCM cluster'
    write(*,*) '    Should be called after setup_quadratic'
    write(*,*) '  anharmonic :'
    write(*,*) '    Collates energies from anharmonic dft runs'
    write(*,*) '    Runs anharmonic calculations'
    write(*,*) '    Should be called after tcm_cluster_run_quadratic'
    write(*,*) '  bs_quadratic :'
    write(*,*) '    Collates bands from anharmonic dft runs'
    write(*,*) '    Runs band structure calculations'
    write(*,*) '    Should be called after tcm_cluster_run_quadratic'
    write(*,*) '  get_kpoints :'
    write(*,*) '    [Help text pending]'
    write(*,*) '  eigenval_vasp_to_bands'
    write(*,*) '    [Help text pending]'
    write(*,*) '  calculate_gap'
    write(*,*) '    [Help text pending]'
  ! test
  elseif (args(1) == 'test') then
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
    temp_real = huge(1.d0)
    temp_real = -1.21313852748395029384e25
    temp_string = temp_real
    write(*,*) temp_real
    write(*,*) char(temp_string)
    temp_real2 = dble(temp_string)
    write(*,*) temp_real2
  ! Wrappers for top-level Fortran
  elseif (args(1) == 'setup_quadratic') then
    call setup_quadratic()
  elseif (args(1) == 'anharmonic') then
    call anharmonic()
  elseif (args(1) == 'bs_quadratic') then
    call bs_quadratic()
  ! Wrappers for subsidiary Fortran 
  elseif (args(1) == 'calculate_gap') then
    call calculate_gap()
  elseif (args(1) == 'combine_forces') then
    call combine_forces(args(2:))
  elseif (args(1) == 'compare_kpoints') then
    call compare_kpoints(args(2:))
  elseif (args(1) == 'construct_finite_displacement') then
    call construct_finite_displacement(args(2:))
  elseif (args(1) == 'construct_matrix_force_cnsts') then
    call construct_matrix_force_cnsts(args(2:))
  elseif (args(1) == 'construct_supercell') then
    call construct_supercell(args(2:))
  elseif (args(1) == 'equilibrium_frac') then
    call equilibrium_frac(args(2:))
  elseif (args(1) == 'fourier_interpolation') then
    call fourier_interpolation(args(2),args(3),args(4),args(5),args(6),args(7),&
      & args(8),args(9),args(10),args(11),args(12))
  elseif (args(1) == 'generate_kgrid') then
    call generate_kgrid(args(2:))
  elseif (args(1) == 'generate_supercells') then
    call generate_supercells(args(2:))
  elseif (args(1) == 'hartree_to_eV') then
    call hartree_to_eV()
  elseif (args(1) == 'lte') then
    tol = dble(args(2))
    tol2 = dble(args(3))
    delta = dble(args(4))
    call lte(tol,tol2,delta,char(args(5)),char(args(6)),char(args(7)), &
      & char(args(8)),char(args(9)),char(args(10)),char(args(11)),     &
      & char(args(12)),char(args(13)),char(args(14)),char(args(15)),   &
      & char(args(16)),char(args(17)),char(args(18)),char(args(19)))
  elseif (args(1) == 'rundft') then
    call rundft(args(2:))
  elseif (args(1) == 'fetch_forces') then
    call fetch_forces(args(2),args(3),int(args(4)),int(args(5)),args(6))
  elseif (args(1) == 'structure_to_dft') then
    if (args(2) == "castep") then
      if (size(args) == 4) then
        call structure_to_dft( dft_code=args(2),              &
                             & structure_sc_filename=args(3), &
                             & output_filename=args(4))
      elseif (size(args) == 5) then
        call structure_to_dft( dft_code=args(2),              &
                             & structure_sc_filename=args(3), &
                             & input_filename=args(4),        &
                             & output_filename=args(5))
      elseif (size(args) == 7) then
        call structure_to_dft( dft_code=args(2),              &
                             & structure_sc_filename=args(3), &
                             & input_filename=args(4),        &
                             & supercell_filename=args(5),    &
                             & path_filename=args(6),         &
                             & output_filename=args(7))
      endif
    elseif (args(2) == "vasp") then
      call structure_to_dft( dft_code=args(2),              &
                           & structure_sc_filename=args(3), &
                           & output_filename=args(4))
    elseif (args(2) == "qe") then
      call structure_to_dft( dft_code=args(2),              &
                           & structure_sc_filename=args(3), &
                           & input_filename=args(4),        &
                           & pseudo_filename=args(5),       &
                           & kpoints_filename=args(6),      &
                           & structure_filename=args(7),    &
                           & output_filename=args(8))
    endif
  elseif (args(1) == 'calculate_symmetry_helper') then
    call calculate_symmetry_helper(args(2:))
  ! wrappers for main shell scripts
  elseif (args(1) == 'setup_harmonic') then
    call system(args(1)//'.sh '//argstring)
  elseif (args(1) == 'tcm_cluster_run_harmonic') then
    call system(args(1)//'.sh '//argstring)
  elseif (args(1) == 'lte_harmonic') then
    call system(args(1)//'.sh '//argstring)
  elseif (args(1) == 'clear_all') then
    call system(args(1)//'.sh '//argstring)
  elseif (args(1) == 'tcm_cluster_run_quadratic') then
    call system(args(1)//'.sh '//argstring)
  elseif (args(1) == 'eigenval_castep_to_bands') then
    call system(args(1)//'.sh '//argstring)
  elseif (args(1) == 'eigenval_vasp_to_bands') then
    call system(args(1)//'.sh '//argstring)
  ! wrappers for subsidiary shell scripts
  elseif (args(1) == 'calculate_symmetry') then
    call system(args(1)//'.sh '//argstring)
  ! wrappers for python scripts
  elseif (args(1) == 'get_kpoints') then
    call system(args(1)//'.py')
  ! unrecognised argument
  else
    write(*,*) char('Unrecognised argument : '//args(1))
  endif
end program

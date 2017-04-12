program caesar
  ! use utility modules
  use string_module
  use io_module
  use utils_module, only : command_line_args, format_path
  
  ! use harmonic modules
  use setup_harmonic_module
  use run_harmonic_module
  use lte_harmonic_module
  
  ! use quadratic modules
  use setup_quadratic_module
  use run_quadratic_module
  use anharmonic_module
  use bs_quadratic_module
  
  ! use testing modules
  use test_copy_harmonic_module
  use test_lte_module
  use test_copy_quadratic_module
  
  ! use misc modules
  use calculate_gap_module
  use hartree_to_eV_module
  
  implicit none
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  type(String)              :: flags_without_arguments
  type(String)              :: flags_with_arguments
  type(CommandLineFlag)     :: flag
  type(String)              :: mode
  
  ! Working directories.
  type(String) :: wd
  type(String) :: cwd
  
  ! Flags.
  type(String) :: input_filename
  logical      :: help
  logical      :: interactive
  
  ! Set terminal width for formatting purposes.
  call update_terminal_width()
  
  ! Get current directory.
  cwd = get_current_directory()
  
  ! Read in command line arguments.
  args = command_line_args()
  if (size(args) < 1) then
    call print_line('No arguments given. For help, call caesar -h')
    stop
  endif
  
  ! Set defaults.
  wd = ''
  input_filename = ''
  help = .false.
  interactive = .false.
  mode = ''
  
  ! Parse command line flags.
  flags_without_arguments = 'hi'
  flags_with_arguments = 'df'
  do
    flag = get_flag(args, flags_without_arguments, flags_with_arguments)
    
    ! No flags remaining.
    if (flag%flag==' ' .and. flag%argument=='') then
      exit
    
    ! An argument which is not a flag.
    elseif (flag%flag==' ') then
      mode = flag%argument
    
    ! Flags without arguments.
    elseif (flag%flag=='h') then
      help = .true.
    
    elseif (flag%flag=='i') then
      interactive = .true.
    
    ! Flags with arguments.
    elseif (flag%flag=='d') then
      wd = flag%argument
    
    elseif (flag%flag=='f') then
      input_filename = flag%argument
    
    endif
  enddo
  
  ! Print help.
  if (help .or. mode == '--help') then
    call print_line('caesar [-h] [-i] [-f input_file] [-d working_directory] &
       &[option]')
    call print_line('')
    call print_line('Flags')
    call print_line('  -h, --help')
    call print_line('      Displays this help text.')
    call print_line('')
    call print_line('  -i')
    call print_line('      Run interactively.')
    call print_line('')
    call print_line('  -f input_file')
    call print_line('      Read settings from specified input file.')
    call print_line('')
    call print_line('  -d working_directory')
    call print_line('      Work in specified directory.')
    call print_line('')
    call print_line('Harmonic options')
    call print_line('  setup_harmonic')
    call print_line('      Sets up harmonic calculation.')
    call print_line('      DFT code choices are: castep.')
    call print_line('  run_harmonic')
    call print_line('      Runs harmonic calculation.')
    call print_line('      Should be called after setup_harmonic.')
    call print_line('  lte_harmonic')
    call print_line('      Runs harmonic calculations.')
    call print_line('      Should be called after run_harmonic.')
    call print_line('')
    call print_line('Quadratic options')
    call print_line('  setup_quadratic')
    call print_line('      Sets up quadratic calculation.')
    call print_line('      DFT code choices are: castep.')
    call print_line('      Should be called after lte_harmonic.')
    call print_line('  run_quadratic')
    call print_line('      Runs quadratic calculation.')
    call print_line('      Should be called after setup_quadratic.')
    call print_line('  anharmonic')
    call print_line('      Runs anharmonic calculations.')
    call print_line('      Should be called after run_quadratic.')
    call print_line('  bs_quadratic')
    call print_line('      Runs band structure calculations.')
    call print_line('      Should be called after run_quadratic.')
    call print_line('')
    call print_line('Utility options')
    call print_line('  hartree_to_eV')
    call print_line('      Provides a Hartree to eV calculator.')
    call print_line('  get_kpoints')
    call print_line('      [Help text pending]')
    call print_line('  calculate_gap')
    call print_line('      [Help text pending]')
    stop
  endif
  
  if (mode=='') then
    call print_line('Error: no mode specified. Call caesar -h for help.')
    stop
  endif
  
  ! Get working directory from user.
  if (wd=='') then
    if (interactive) then
      call print_line('')
      call print_line('Where is the working directory?')
      wd = format_path(read_line_from_user(), cwd)
    else
      wd = '.'
    endif
  endif
  
  ! Run main program, depending on mode.
  if (mode == 'test') then
    call system_call(str('mkdir this_is_another_dir'))
  ! Wrappers for top-level Fortran.
  elseif (mode == 'setup_harmonic') then
    call setup_harmonic(wd)
  elseif (mode == 'run_harmonic') then
    call run_harmonic(wd,cwd)
  elseif (mode == 'lte_harmonic') then
    call lte_harmonic(wd)
  elseif (mode == 'setup_quadratic') then
    call setup_quadratic(wd,cwd)
  elseif (mode == 'run_quadratic') then
    call run_quadratic(wd,cwd)
  elseif (mode == 'anharmonic') then
    call anharmonic(wd)
  elseif (mode == 'bs_quadratic') then
    call bs_quadratic(wd)
  ! Wrappers for subsidiary Fortran 
  elseif (mode == 'calculate_gap') then
    call calculate_gap()
  elseif (mode == 'hartree_to_eV') then
    call hartree_to_eV()
  ! Wrappers for testing modules
  elseif (mode == 'test_copy_harmonic') then
    call test_copy_harmonic(wd,cwd)
  elseif (mode == 'test_lte') then
    call test_lte(wd,cwd)
  elseif (mode == 'test_copy_quadratic') then
    call test_copy_quadratic(wd,cwd)
  ! unrecognised argument
  else
    call print_line('Error: Unrecognised argument: '//mode)
  endif
end program

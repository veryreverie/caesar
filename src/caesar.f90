! ======================================================================
! The main program of Caesar.
!
! Processes user inputs and calls a subsidiary function.
! Call Caesar -h for details.
! ======================================================================
program caesar
  ! Use utility modules.
  use string_module
  use io_module
  use dictionary_module
  use utils_module, only : command_line_args
  use help_module
  use process_arguments_module
  
  ! Use harmonic modules.
  use setup_harmonic_module
  use run_harmonic_module
  use calculate_harmonic_module
  
  ! Use quadratic modules.
  use setup_quadratic_module
  use run_quadratic_module
  use anharmonic_module
  use bs_quadratic_module
  
  ! Use anharmonic modules.
  use setup_anharmonic_module
  use run_anharmonic_module
  use calculate_anharmonic_module
  
  ! Use testing modules.
  use test_module
  use linear_algebra_test_module
  use setup_harmonic_test_module
  use calculate_harmonic_test_module
  use test_copy_quadratic_module
  
  ! Use misc modules.
  use calculate_gap_module
  use hartree_to_eV_module
  
  implicit none
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  
  ! The first command line argument, the mode.
  type(String)              :: mode
  
  ! The chosen subroutine, and the keywords it accepts.
  procedure(MainSubroutine), pointer     :: main_subroutine => null ()
  type(KeywordData),         allocatable :: keywords(:)
  
  ! Processed input arguments.
  type(Dictionary)          :: arguments
  
  ! The working directory.
  type(String) :: wd
  
  ! Temporary variables.
  type(String)              :: filename
  
  ! --------------------------------------------------
  ! Set IO variables for formatting and file parsing purposes.
  ! --------------------------------------------------
  call set_global_io_variables()
  
  ! --------------------------------------------------
  ! Read in command line arguments and process the mode.
  ! --------------------------------------------------
  
  ! Read in arguments.
  args = command_line_args()
  
  ! Error: no arguments given.
  if (size(args) < 2) then
    call print_line('No arguments given. For help, call caesar -h')
    stop
  endif
  
  ! Read in mode.
  mode = lower_case(args(2))
  
  ! Error: mode only contains one character.
  if (len(mode) < 2) then
    call print_line('Error: unrecognised mode: '//mode)
    call print_line('Call caesar -h for help.')
    stop
  
  ! Help calls, both correct and malformed.
  elseif (mode == '-h' .or. mode == '--help' .or. mode == 'help') then
    if (size(args) == 2) then
      call help()
    else
      call print_line('For keyword-specific help, please also specify a mode, &
         &e.g. "caesar setup_harmonic -h dft_code"')
    endif
    stop
  elseif (slice(mode,1,2) == '-h') then
    call print_line('For keyword-specific help, please also specify a mode, &
       &e.g. "caesar setup_harmonic -h dft_code"')
    stop
  
  ! Normal inputs. Fetch keywords and set subprocess.
  elseif (mode == 'test') then
    keywords = test_keywords()
    main_subroutine => test
  elseif (mode == 'setup_harmonic') then
    keywords = setup_harmonic_keywords()
    main_subroutine => setup_harmonic
  elseif (mode == 'run_harmonic') then
    keywords = run_harmonic_keywords()
    main_subroutine => run_harmonic
  elseif (mode == 'calculate_harmonic') then
    keywords = calculate_harmonic_keywords()
    main_subroutine => calculate_harmonic
  elseif (mode == 'setup_quadratic') then
    keywords = setup_quadratic_keywords()
    main_subroutine => setup_quadratic
  elseif (mode == 'run_quadratic') then
    keywords = run_quadratic_keywords()
    main_subroutine => run_quadratic
  elseif (mode == 'anharmonic') then
    keywords = anharmonic_keywords()
    main_subroutine => anharmonic
  elseif (mode == 'bs_quadratic') then
    keywords = bs_quadratic_keywords()
    main_subroutine => bs_quadratic
  elseif (mode == 'setup_anharmonic') then
    keywords = setup_anharmonic_keywords()
    main_subroutine => setup_anharmonic
  elseif (mode == 'run_anharmonic') then
    keywords = run_anharmonic_keywords()
    main_subroutine => run_anharmonic
  elseif (mode == 'calculate_anharmonic') then
    keywords = calculate_anharmonic_keywords()
    main_subroutine => calculate_anharmonic
  elseif (mode == 'linear_algebra_test') then
    keywords = linear_algebra_test_keywords()
    main_subroutine => linear_algebra_test
  elseif (mode == 'setup_harmonic_test') then
    keywords = setup_harmonic_test_keywords()
    main_subroutine => setup_harmonic_test
  elseif (mode == 'calculate_harmonic_test') then
    keywords = calculate_harmonic_test_keywords()
    main_subroutine => calculate_harmonic_test
  elseif (mode == 'test_copy_quadratic') then
    keywords = test_copy_quadratic_keywords()
    main_subroutine => test_copy_quadratic
  elseif (mode == 'calculate_gap') then
    keywords = calculate_gap_keywords()
    main_subroutine => calculate_gap
  elseif (mode == 'hartree_to_ev') then
    keywords = hartree_to_eV_keywords()
    main_subroutine => hartree_to_eV
  
  ! Erroneous inputs.
  elseif (slice(mode,1,1) == '-') then
    call print_line('Error: The first argument should be the mode, and not a &
       &flag or keyword.')
    call print_line('Call caesar -h for help.')
    stop
  else
    call print_line('Error: unrecognised mode: '//mode)
    call print_line('Call caesar -h for help.')
    stop
  endif
  
  ! --------------------------------------------------
  ! Process arguments and flags.
  ! --------------------------------------------------
  arguments = process_arguments(args,keywords)
  
  ! --------------------------------------------------
  ! Handle help calls.
  ! --------------------------------------------------
  if (arguments%is_set('help')) then
    call help(arguments%value('help'), mode, keywords)
    stop
  endif
  
  ! --------------------------------------------------
  ! Write settings to file.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  filename = wd//'/'//mode//'.used_settings'
  call arguments%write_file(filename)
  call print_line('')
  call print_line('Settings written to file '//filename)
  
  ! --------------------------------------------------
  ! Run main program with input arguments, depending on mode.
  ! --------------------------------------------------
  call main_subroutine(arguments)
end program

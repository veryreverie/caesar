! ======================================================================
! The main program of Caesar.
!
! Processes user inputs and calls a subsidiary function.
! Call Caesar -h for details.
! ======================================================================
program caesar
  ! Use utility modules.
  use common_module
  
  ! Use convergence modules.
  use converge_cutoff_and_kpoints_module
  use plot_cutoff_and_kpoints_module
  
  ! Use harmonic modules.
  use setup_harmonic_module
  use run_harmonic_module
  use calculate_normal_modes_module
  use plot_normal_modes_module
  use calculate_harmonic_observables_module
  use plot_dos_and_dispersion_module
  use plot_thermodynamic_variables_module
  
  ! Use anharmonic modules.
  use setup_anharmonic_module
  use run_anharmonic_module
  use calculate_potential_module
  use calculate_states_module
  
  ! Use testing modules.
  use test_module
  use linear_algebra_test_module
  
  ! Use misc modules.
  use hartree_to_eV_module
  use version_module
  implicit none
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  
  ! The first command line argument, the mode.
  type(String) :: mode
  
  ! The chosen subroutine, and the keywords it accepts.
  type(CaesarModes)                      :: caesar_modes
  type(CaesarMode)                       :: caesar_mode
  procedure(MainSubroutine), pointer     :: main_subroutine => null ()
  type(KeywordData),         allocatable :: keywords(:)
  
  ! Processed input arguments.
  type(Dictionary) :: arguments
  
  ! The working directory.
  type(String) :: wd
  
  ! The output file where stdout is piped if -o is specified.
  type(OFile) :: output_file
  
  ! Temporary variables.
  type(String)              :: filename
  
  ! --------------------------------------------------
  ! Set IO variables for formatting and file parsing purposes.
  ! --------------------------------------------------
  call set_io_settings()
  
  ! --------------------------------------------------
  ! Read in mode interfaces.
  ! --------------------------------------------------
  ! Normal inputs. Fetch keywords and set subprocess.
  caesar_modes = CaesarModes([           &
     & hartree_to_ev(),                  &
     & test(),                           &
     & converge_cutoff_and_kpoints(),    &
     & plot_cutoff_and_kpoints(),        &
     & setup_harmonic(),                 &
     & run_harmonic(),                   &
     & calculate_normal_modes(),         &
     & plot_normal_modes(),              &
     & calculate_harmonic_observables(), &
     & plot_dos_and_dispersion(),        &
     & plot_thermodynamic_variables(),   &
     & setup_anharmonic(),               &
     & run_anharmonic(),                 &
     & calculate_potential(),            &
     & calculate_states()                &
     &])
  
  ! --------------------------------------------------
  ! Read in command line arguments and process the mode.
  ! --------------------------------------------------
  
  ! Read in arguments.
  args = command_line_args()
  
  ! Error: no arguments given.
  if (size(args) < 2) then
    call print_line(colour('Error: no arguments given.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop
  endif
  
  ! Read in mode.
  mode = lower_case(args(2))
  
  ! Error: mode only contains one character.
  if (len(mode) < 2) then
    call print_line(colour('Error: unrecognised mode: ','red')//mode)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop
  
  ! Help calls, both correct and malformed.
  elseif (mode == '-h' .or. mode == '--help' .or. mode == 'help') then
    if (size(args) == 2) then
      call help(caesar_modes)
    else
      call print_line(colour('Error: no mode specified.','red'))
      call print_line('For keyword-specific help, please also specify a mode, &
         &as')
      call print_line(colour('caesar [mode] -h [keyword]','white')//'.')
    endif
    stop
  elseif (slice(mode,1,2) == '-h') then
    call print_line(colour('Error: no mode specified.','red'))
    call print_line('For keyword-specific help, please also specify a mode, &
       &as')
    call print_line(colour('caesar [mode] -h [keyword]','white')//'.')
    stop
  
  ! Version calls.
  elseif (mode=='--version') then
    call print_version()
    stop
  
  ! Erroneous inputs.
  elseif (slice(mode,1,1) == '-') then
    call print_line(colour('Error: The first argument should be the mode, &
       &and not a flag or keyword.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop
  
  ! Normal inputs.
  else
    caesar_mode = caesar_modes%mode(mode)
    keywords = caesar_mode%keywords
    main_subroutine => caesar_mode%main_subroutine
  endif
  
  ! --------------------------------------------------
  ! Process arguments and flags.
  ! --------------------------------------------------
  arguments = process_arguments(args,keywords)
  
  ! --------------------------------------------------
  ! Set output file, if appropriate.
  ! --------------------------------------------------
  if (arguments%is_set('output_file')) then
    output_file = OFile(arguments%value('output_file'))
    call output_file%make_stdout()
  endif
  
  ! --------------------------------------------------
  ! Handle help calls.
  ! --------------------------------------------------
  if (arguments%is_set('help')) then
    call help(arguments%value('help'), mode, caesar_modes)
    stop
  endif
  
  ! --------------------------------------------------
  ! Print version and quit, if requested.
  ! --------------------------------------------------
  if (arguments%is_set('version')) then
    call print_version()
    stop
  endif
  
  ! --------------------------------------------------
  ! Write settings to file.
  ! --------------------------------------------------
  if (.not. caesar_mode%suppress_settings_file) then
    wd = arguments%value('working_directory')
    filename = wd//'/'//mode//'.used_settings'
    call arguments%write_file(filename)
  endif
  
  ! --------------------------------------------------
  ! Run main program with input arguments, depending on mode.
  ! --------------------------------------------------
  call main_subroutine(arguments)
end program

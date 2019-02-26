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
  use harmonic_module
  
  ! Use anharmonic modules.
  use anharmonic_module
  
  ! Use Castep modules.
  use converge_harmonic_frequencies_module
  use plot_harmonic_convergence_module
  
  ! Use testing modules.
  use check_counter_module
  use test_module
  use linear_algebra_test_module
  use update_basis_functions_module
  
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
  
  ! The output file where stdout is piped if -o is specified.
  type(OFile) :: output_file
  
  ! Temporary variables.
  type(String) :: filename
  
  ! --------------------------------------------------
  ! Call startup subroutines,
  !    which must be called before anything else happens.
  ! --------------------------------------------------
  call startup_utils()
  call startup_anharmonic()
  
  ! --------------------------------------------------
  ! Read in mode interfaces.
  ! --------------------------------------------------
  ! Normal inputs. Fetch keywords and set subprocess.
  caesar_modes = CaesarModes([             &
     & test(),                             &
     & check_counter(),                    &
     & hartree_to_ev(),                    &
     & converge_cutoff_and_kpoints(),      &
     & plot_cutoff_and_kpoints(),          &
     & setup_harmonic(),                   &
     & run_harmonic(),                     &
     & calculate_normal_modes(),           &
     & plot_normal_modes(),                &
     & calculate_harmonic_observables(),   &
     & plot_dos_and_dispersion(),          &
     & plot_thermodynamic_variables(),     &
     & converge_qpoint_grid(),             &
     & setup_anharmonic(),                 &
     & run_anharmonic(),                   &
     & calculate_potential(),              &
     & map_anharmonic_modes(),             &
     & plot_anharmonic_modes(),            &
     & map_potential(),                    &
     & plot_potential_map(),               &
     & map_vscf_modes(),                   &
     & plot_vscf_modes(),                  &
     & calculate_anharmonic_observables(), &
     & converge_harmonic_frequencies(),    &
     & plot_harmonic_convergence(),        &
     & plot_vscf_states(),                 &
     & update_basis_functions()            ])
  
  ! --------------------------------------------------
  ! Read in command line arguments and process the mode.
  ! --------------------------------------------------
  
  ! Read in arguments.
  args = command_line_args()
  
  ! Error: no arguments given.
  if (size(args) < 2) then
    call print_line(colour('Error: no arguments given.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop 1
  endif
  
  ! Read in mode.
  mode = lower_case(args(2))
  
  ! Error: mode only contains one character.
  if (len(mode) < 2) then
    call print_line(colour('Error: unrecognised mode: ','red')//mode)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop 1
  
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
    stop 1
  elseif (slice(mode,1,2) == '-h') then
    call print_line(colour('Error: no mode specified.','red'))
    call print_line('For keyword-specific help, please also specify a mode, &
       &as')
    call print_line(colour('caesar [mode] -h [keyword]','white')//'.')
    stop 1
  
  ! Version calls.
  elseif (mode=='--version') then
    call print_version()
    stop 1
  
  ! Erroneous inputs.
  elseif (slice(mode,1,1) == '-') then
    call print_line(colour('Error: The first argument should be the mode, &
       &and not a flag or keyword.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    stop 1
  
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
    stop 1
  endif
  
  ! --------------------------------------------------
  ! Print version and quit, if requested.
  ! --------------------------------------------------
  if (arguments%is_set('version')) then
    call print_version()
    stop 1
  endif
  
  ! --------------------------------------------------
  ! Set working directory.
  ! --------------------------------------------------
  call set_working_directory(arguments%value('working_directory'))
  
  ! --------------------------------------------------
  ! Write settings to file.
  ! --------------------------------------------------
  if (.not. caesar_mode%suppress_settings_file) then
    filename = mode//'.used_settings'
    call arguments%write_file(filename)
  endif
  
  ! --------------------------------------------------
  ! Run main program with input arguments, depending on mode.
  ! --------------------------------------------------
  call main_subroutine(arguments)
end program

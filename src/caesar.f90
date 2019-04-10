! ======================================================================
! The main program of Caesar.
!
! Processes user inputs and calls a subsidiary function.
! Call Caesar -h for details.
! ======================================================================
program caesar
  use common_module
  use harmonic_module
  use anharmonic_module
  use dft_module
  use testing_module
  use version_module
  implicit none
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  
  ! The first command line argument, the mode.
  type(String) :: mode
  
  ! The chosen subroutine, and the keywords it accepts.
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
  call startup_harmonic()
  call startup_anharmonic()
  call startup_dft()
  call startup_testing()
  
  ! --------------------------------------------------
  ! Read in command line arguments and process the mode.
  ! --------------------------------------------------
  
  ! Read in arguments.
  args = command_line_args()
  
  ! Error: no arguments given.
  if (size(args) < 2) then
    call print_line(colour('Error: no arguments given.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  endif
  
  ! Read in mode.
  mode = lower_case(args(2))
  
  ! Error: mode only contains one character.
  if (len(mode) < 2) then
    call print_line(colour('Error: unrecognised mode: ','red')//mode)
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  
  ! Help calls, both correct and malformed.
  elseif (mode == '-h' .or. mode == '--help' .or. mode == 'help') then
    if (size(args) == 2) then
      call help()
    else
      call print_line(colour('Error: no mode specified.','red'))
      call print_line('For keyword-specific help, please also specify a mode, &
         &as')
      call print_line(colour('caesar [mode] -h [keyword]','white')//'.')
    endif
    call quit()
  elseif (slice(mode,1,2) == '-h') then
    call print_line(colour('Error: no mode specified.','red'))
    call print_line('For keyword-specific help, please also specify a mode, &
       &as')
    call print_line(colour('caesar [mode] -h [keyword]','white')//'.')
    call quit()
  
  ! Version calls.
  elseif (mode=='--version') then
    call print_version()
    call quit()
  
  ! Erroneous inputs.
  elseif (slice(mode,1,1) == '-') then
    call print_line(colour('Error: The first argument should be the mode, &
       &and not a flag or keyword.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  
  ! Normal inputs.
  else
    caesar_mode = CaesarMode(mode)
    keywords = caesar_mode%keywords
    main_subroutine => caesar_mode%main_subroutine
  endif
  
  ! --------------------------------------------------
  ! Process arguments and flags.
  ! --------------------------------------------------
  arguments = Dictionary(args,keywords)
  
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
    call help(arguments%value('help'), mode)
    call quit()
  endif
  
  ! --------------------------------------------------
  ! Print version and quit, if requested.
  ! --------------------------------------------------
  if (arguments%is_set('version')) then
    call print_version()
    call quit()
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

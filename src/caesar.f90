! Caesar; a utility for calculating the vibrational free energy of periodic crystals.
! Copyright (C) 2021 Mark Johnson
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!> The `caesar` executable. Processes user inputs and calls the relevant
!>    [[ProgramMode(type)]].
!> Call `caesar -h` for details.
program caesar
  use caesar_common_module
  use caesar_harmonic_module
  use caesar_anharmonic_module
  use caesar_dft_module
  use caesar_testing_module
  use caesar_version_module
  implicit none
  
  ! Command line arguments.
  type(String), allocatable :: args(:)
  
  ! The first command line argument, the mode.
  type(String) :: mode
  
  ! The chosen subroutine, and the keywords it accepts.
  type(ProgramModes)                     :: program_modes
  type(ProgramMode)                      :: program_mode
  procedure(MainSubroutine), pointer     :: main_subroutine => null ()
  type(KeywordData),         allocatable :: keywords(:)
  
  ! Processed input arguments.
  type(Dictionary) :: arguments
  
  ! The output file where stdout is piped if -o is specified.
  type(OFile) :: output_file
  
  ! Temporary variables.
  type(String) :: filename
  
  ! --------------------------------------------------
  ! Collect together the program modes.
  ! --------------------------------------------------
  program_modes = ProgramModes([ harmonic_modes(),   &
                               & anharmonic_modes(), &
                               & dft_modes(),        &
                               & testing_modes()     ])
  
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
      call help(program_modes)
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
  
  ! Copyright calls.
  elseif (mode=='--copyright') then
    call print_copyright()
    call quit()
  
  ! Erroneous inputs.
  elseif (slice(mode,1,1) == '-') then
    call print_line(colour('Error: The first argument should be the mode, &
       &and not a flag or keyword.','red'))
    call print_line('Call '//colour('caesar -h','white')//' for help.')
    call quit()
  
  ! Normal inputs.
  else
    program_mode = program_modes%mode(mode)
    keywords = program_mode%keywords
    main_subroutine => program_mode%main_subroutine
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
    call help(arguments%value('help'), program_mode)
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
  call set_current_working_directory(arguments%value('working_directory'))
  
  ! --------------------------------------------------
  ! Write settings to file.
  ! --------------------------------------------------
  if (.not. program_mode%suppress_settings_file) then
    filename = mode//'.used_settings'
    call arguments%write_file(filename)
  endif
  
  ! --------------------------------------------------
  ! Run main program with input arguments, depending on mode.
  ! --------------------------------------------------
  call main_subroutine(arguments)
end program

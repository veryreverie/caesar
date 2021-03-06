submodule (caesar_help_module) caesar_help_submodule
  use caesar_arguments_module
contains

module procedure help_default
  call print_line('Caesar; a utility for calculating the vibrational free &
     &energy of periodic crystals.')
  call print_line('Copyright (C) 2021 Mark Johnson')
  call print_line('This program comes with ABSOLUTELY NO WARRANTY. &
     &This is free software, and you are welcome to redistribute it &
     &under certain conditions; type "caesar --copyright" for details.')
  call print_line('')
  call print_line('Usage:')
  call print_line('')
  call print_line( colour('caesar','white')                          // &
                 & colour(' mode ','cyan')                           // &
                 & '['//colour('-h','white')//' [keyword]] '         // &
                 & '['//colour('-i','white')//'] '                   // &
                 & '['//colour('-f','white')//' input_file] '        // &
                 & '['//colour('-d','white')//' working_directory] ' // &
                 & '['//colour('-o','white')//' output_file] '       // &
                 & '[--options]')
  call print_line('')
  call print_line('')
  call print_line('Accepted arguments to all modes:')
  call print_line('')
  call print_line( colour('-h','white')//' [keyword] | '// &
                 & colour('--help','white')//' [keyword]',   &
                 & settings=PrintSettings(indent=2))
  call print_line( colour('caesar -h','white')//         &
                 & '              : Displays this help text.', &
                 & settings=PrintSettings(indent=4))
  call print_line( colour('caesar','white')//colour(' mode','cyan')// &
                 & colour(' -h','white')//'         : Displays help text &
                 &relevant to the specified '//colour('mode','cyan')//'.', &
                 & settings=PrintSettings(indent=4))
  call print_line( colour('caesar','white')//colour(' mode','cyan')// &
                 & colour(' -h','white')//' keyword : displays help text &
                 &relevant to the specified keyword, in the context of the &
                 &specified '//colour('mode','cyan')//'.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('-i','white')//' | '// &
                 & colour('--interactive','white'), &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Runs interactively, prompting the user to review and &
                 &set all options.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('-f','white')//' filename | '// &
                 & colour('--input_file','white')//' filename', &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Reads additional settings from specified file.', &
                 & settings=PrintSettings(indent=4))
  call print_line( 'These should be of the form:', &
                 & settings=PrintSettings(indent=4))
  call print_line( 'keyword1 argument', &
                 & settings=PrintSettings(indent=6))
  call print_line( 'keyword2 argument1 argument2 argument3', &
                 & settings=PrintSettings(indent=6))
  call print_line( 'Keywords are the same as command-line --keywords.', &
                 & settings=PrintSettings(indent=4))
  call print_line( 'The "--" prefix. should not be given.', &
                 & settings=PrintSettings(indent=4))
  call print_line( 'The keywords "filename", "interactive" and "help" &
                 &should not be specified in a file.',               &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('-d','white')//' directory_name | '// &
                 & colour('--working_directory','white')//' directory_name', &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Specifies the directory where Caesar should work.', &
                 & settings=PrintSettings(indent=4))
  call print_line( 'All files and folders will be created here.', &
                 & settings=PrintSettings(indent=4))
  call print_line( 'This is also where any run scripts will be called.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('-o','white')//' filename | '// &
                 & colour('--output_file','white')//' filename', &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Specifies a file to which all terminal output will be &
                 &written. This also disables terminal formatting, so should &
                 &be favoured over piping to file. If unset, terminal output &
                 &will go to the terminal.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('--random_seed','white')//' seed', &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Specifies the seed which will be used to initialise any &
                 &random number generation, allowing computations to be &
                 &repeated exactly. If unset, the seed will be set to the &
                 &current time in milliseconds.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('--version','white'), &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Displays version information and quits.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line( colour('--copyright','white'), &
                 & settings=PrintSettings(indent=2))
  call print_line( 'Displays copyright information and quits.', &
                 & settings=PrintSettings(indent=4))
  call print_line('')
  call print_line('')
  call print_line('Accepted '//colour('modes','cyan')//':')
  
  call program_modes%print_help()
  
  call print_line('')
  call print_line('Return code:')
  call print_line('')
  call print_line( 'Caesar will return 0 on successful termination, and 1 &
                 &otherwise.',                                            &
                 & settings=PrintSettings(indent=2))
  call print_line('')
  call print_line('Suggested usage:')
  call print_line('')
  call print_line( 'For harmonic phonon calculations:', &
                 & settings=PrintSettings(indent=2)     )
  call print_line( colour('caesar','white')//' '//  &
                 & colour('setup_harmonic','cyan'), &
                 & settings=PrintSettings(indent=4) )
  call print_line( colour('caesar','white')//' '//  &
                 & colour('run_harmonic','cyan'),   &
                 & settings=PrintSettings(indent=4) )
  call print_line( colour('caesar','white')//' '//          &
                 & colour('calculate_normal_modes','cyan'), &
                 & settings=PrintSettings(indent=4)         )
  call print_line( colour('caesar','white')//' '//                  &
                 & colour('calculate_harmonic_observables','cyan'), &
                 & settings=PrintSettings(indent=4)                 )
  call print_line('')
  call print_line( 'For anharmonic phonon calculations (after harmonic &
                 &calculations):',                                     &
                 & settings=PrintSettings(indent=2)                    )
  call print_line( colour('caesar','white')//' '//    &
                 & colour('setup_anharmonic','cyan'), &
                 & settings=PrintSettings(indent=4)   )
  call print_line( colour('caesar','white')//' '//  &
                 & colour('run_anharmonic','cyan'), &
                 & settings=PrintSettings(indent=4) )
  call print_line( colour('caesar','white')//' '//       &
                 & colour('calculate_potential','cyan'), &
                 & settings=PrintSettings(indent=4)      )
  call print_line( colour('caesar','white')//' '//                    &
                 & colour('calculate_anharmonic_observables','cyan'), &
                 & settings=PrintSettings(indent=4)                   )
end procedure

module procedure help_specific
  type(KeywordData), allocatable :: keywords(:)
  
  integer :: i
  
  keywords = program_mode%keywords
  
  if (keyword=='') then
    call program_mode%print_help()
    call print_line('')
    call print_line('')
    call print_line('Accepted keywords:')
    do i=1,size(keywords)
      call keywords(i)%print_help()
    enddo
    call quit()
  else
    do i=1,size(keywords)
      if (keywords(i)%keyword()==keyword) then
        call keywords(i)%print_help()
        call quit()
      endif
    enddo
    
    call print_line('')
    call print_line('Keyword '//keyword//' not recognised. For a list of &
       &keywords associated with mode '//program_mode%mode_name//', call:')
    call print_line('  caesar '//program_mode%mode_name//' -h')
  endif
end procedure

module procedure print_copyright
  call print_line('Caesar; a utility for calculating the vibrational free &
     &energy of periodic crystals.')
  call print_line('Copyright (C) 2021 Mark Johnson')
  call print_line('')
  call print_line('This program is free software: you can redistribute it &
     &and/or modify it under the terms of the GNU Lesser General Public &
     &License as published by the Free Software Foundation, either version 3 &
     &of the License, or (at your option) any later version.')
  call print_line('')
  call print_line('This program is distributed in the hope that it will be &
     &useful, but WITHOUT ANY WARRANTY; without even the implied warranty of &
     &MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU &
     &Lesser General Public License for more details.')
  call print_line('')
  call print_line('You should have received a copy of the GNU Lesser General &
     &Public License along with this program.  If not, see &
     &<https://www.gnu.org/licenses/>.')
end procedure
end submodule

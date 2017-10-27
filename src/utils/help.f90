! ======================================================================
! Handles input keywords and helptext.
! ======================================================================
module help_module
  use constants_module, only : dp
  use string_module
  use io_module
  use keyword_module
  implicit none
  
  private
  
  ! Prints help text.
  public :: help
  
  interface help
    module procedure help_default
    module procedure help_specific
  end interface
  
contains

! ----------------------------------------------------------------------
! Prints helptext.
! ----------------------------------------------------------------------
  
! Print default helptext.
subroutine help_default(caesar_modes)
  use caesar_modes_module
  implicit none
  
  type(CaesarModes), intent(in) :: caesar_modes
  
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
  call print_line('Accepted Flags:')
  call print_line('')
  call print_line( colour('  -h','white')//' [keyword] | '// &
                 & colour('--help','white')//' [keyword]')
  call print_line( colour('      caesar -h','white')// &
                 & '              : Displays this help text.')
  call print_line( colour('      caesar','white')//colour(' mode','cyan')// &
                 & colour(' -h','white')//'         : Displays help text &
                 &relevant to the specified '//colour('mode','cyan')//'.')
  call print_line( colour('      caesar','white')//colour(' mode','cyan')// &
                 & colour(' -h','white')//' keyword : displays help text &
                 &relevant to the specified keyword, in the context of the &
                 &specified '//colour('mode','cyan')//'.')
  call print_line('')
  call print_line( colour('  -i','white')//' | '// &
                 & colour('--interactive','white'))
  call print_line('      Runs interactively, prompting the user to review and &
     &set all options.')
  call print_line('')
  call print_line( colour('  -f','white')//' filename | '// &
                 & colour('--input_file','white')//' filename')
  call print_line('      Reads additional settings from specified file.')
  call print_line('      These should be of the form:')
  call print_line('         keyword1 argument')
  call print_line('         keyword2 argument1 argument2 argument3')
  call print_line('      Keywords are the same as command-line --keywords.')
  call print_line('      The "--" prefix. should not be given.')
  call print_line('      The keywords "filename", "interactive" and "help" &
     &should not be specified in a file.')
  call print_line('')
  call print_line( colour('  -d','white')//' directory_name | '// &
                 & colour('--working_directory','white')//' directory_name')
  call print_line('      Specifies the directory where Caesar should work.')
  call print_line('      All files and folders will be created here.')
  call print_line('      This is also where any run scripts will be called.')
  call print_line('')
  call print_line( colour('  -o','white')//' filename | '// &
                 & colour('--output_file','white')//' filenamename')
  call print_line('      Specifies a file to which all terminal output will &
     &be written. This also disables terminal formatting, so should be &
     &favoured over piping to file. If unset, terminal output will go to the &
     &terminal.')
  call print_line('')
  call print_line('')
  call print_line('Accepted '//colour('modes','cyan')//':')
  
  call caesar_modes%print_help()
end subroutine

! Prints the helptext for a particular mode or keyword.
subroutine help_specific(keyword,mode,caesar_modes)
  use caesar_modes_module
  implicit none
  
  type(String),      intent(in) :: keyword
  type(String),      intent(in) :: mode
  type(CaesarModes), intent(in) :: caesar_modes
  
  type(CaesarMode)               :: caesar_mode
  type(KeywordData), allocatable :: keywords(:)
  
  integer :: i
  
  caesar_mode = caesar_modes%mode(mode)
  keywords = caesar_mode%keywords
  
  if (keyword=='') then
    call caesar_mode%print_help()
    call print_line('')
    call print_line('')
    call print_line('Accepted keywords:')
    do i=1,size(keywords)
      call keywords(i)%print_help()
    enddo
    stop
  else
    do i=1,size(keywords)
      if (keywords(i)%keyword==keyword) then
        call keywords(i)%print_help()
        stop
      endif
    enddo
    
    call print_line('')
    call print_line('Keyword '//keyword//' not recognised. For a list of &
       &keywords associated with mode '//mode//', call:')
    call print_line('  caesar '//mode//' -h')
  endif
end subroutine
end module

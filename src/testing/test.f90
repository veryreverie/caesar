! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function test_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = test_keywords()
  output%main_subroutine => test
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  use dictionary_module
  use ifile_module
  use ofile_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(String) :: ten
  type(String) :: line
  
  wd = arguments%value('working_directory')
  
  ten = '1234567890'
  
  line = ten//ten//ten//ten//ten//ten//ten//ten//'1234'
  
  call print_line(colour('black','black'))
  call print_line(colour('red','red'))
  call print_line(colour('green','green'))
  call print_line(colour('yellow','yellow'))
  call print_line(colour('blue','blue'))
  call print_line(colour('magenta','magenta'))
  call print_line(colour('cyan','cyan'))
  call print_line(colour('light gray','light gray'))
  call print_line(colour('dark gray','dark gray'))
  call print_line(colour('light red','light red'))
  call print_line(colour('light green','light green'))
  call print_line(colour('light yellow','light yellow'))
  call print_line(colour('light blue','light blue'))
  call print_line(colour('light magenta','light magenta'))
  call print_line(colour('light cyan','light cyan'))
  call print_line(colour('white','white'))
  
  call print_line('This is a really really really long long long sentence, &
     &which seems like it should have stopped some time ago.')
  
  call print_line('This is a '//colour('really','red')//' '// &
     &colour('really','red')//' '//colour('really','red')//' '// &
     &colour('long','blue')//' '//colour('long','blue')//' '// &
     &colour('long','blue')//' sentence, &
     &which seems like it should have stopped some time ago.')
  
  call print_line(line//'overflow')
  call print_line(colour(line,'red')//'overflow')
  call print_line(colour(line,'red')//colour('overflow','blue'))
  call print_line(slice(line,1,70)//' '//line)
  call print_line(slice(line,1,70)//' '//line//'overflow '// &
     & colour('sameline','green')//' '//line)
  call print_line('hi                                               &
  &                                   there')
  call print_line( 'Specifies a file to which all terminal output will be &
                 &written. This also disables terminal formatting, so should &
                 &be favoured over piping to file. If unset, terminal output &
                 &will go to the terminal.', &
                 & indent=4)
  
end subroutine
end module

! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test
  
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'test'
  output%description = 'Runs temporary code for testing purposes.'
  output%keywords = [KeywordData::]
  output%main_subroutine => test_subroutine
  output%suppress_from_helptext = .true.
end function

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: thing
  
  thing = 'abcoisvsniov'
  
  call print_line(thing)
  call print_line(replace(thing,'S','a'))
  call print_line(thing)
end subroutine
end module

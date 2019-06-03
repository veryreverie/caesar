! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  use anharmonic_module
  
  implicit none
  
  private
  
  public :: startup_test
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_test()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'test'
  mode%description = 'Runs temporary code for testing purposes.'
  mode%keywords = [KeywordData::]
  mode%main_subroutine => test_subroutine
  mode%suppress_from_helptext = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: a
  
  a = 'this is a string'
  call print_line(token(a,2))
  call print_line(token(a,1,delimiter='r'))
  call print_lines(tokens(a,2,4))
  call print_lines(tokens(a,3))
  call print_lines(tokens(a,last=3))
  call print_lines(tokens(a,4,18,'x'))
end subroutine
end module

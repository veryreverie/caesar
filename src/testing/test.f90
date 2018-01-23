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
  
  use fraction_module
  use fraction_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(IntFraction) :: a
  type(IntFraction) :: b
  type(IntFraction) :: c
  type(IntFraction) :: d
  
  wd = arguments%value('working_directory')
  
  a = 3
  a = a / 6
  call print_line(a)
  
end subroutine
end module

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

subroutine test(arguments)
  use dictionary_module
  use linear_algebra_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: a
  
  a = 'f'
  call print_line(lgcl(a))
  
  a = 't'
  call print_line(lgcl(a))
  
  a = 'F'
  call print_line(lgcl(a))
  
  a = 'T'
  call print_line(lgcl(a))
  
  a = 'False'
  call print_line(lgcl(a))
  
  a = 'TrUe'
  call print_line(lgcl(a))
  
  a = 'hi'
  call print_line(lgcl(a))
  
  call print_line('END TEST')
end subroutine
end module

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

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  use dictionary_module
  use stringable_example_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(StringableExample) :: a
  type(String)            :: b
  
  a%contents = 1
  
  b = a
  
  call print_line(b)
  call print_line('a: '//a)
  call print_line(a//': a')
  call print_line(1//': '//a)
  call print_line(a//': '//1)
  
end subroutine
end module

! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
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
  use single_mode_displacement_module
  use mode_displacement_module
  use univariate_module
  use monomial_module
  use polynomial_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  integer, allocatable :: list1(:)
  integer, allocatable :: list2(:)
  integer, allocatable :: list3(:)
  
  integer :: i
  
  wd = arguments%value('working_directory')
  
  list1 = [1,2,3]
  list2 = [integer::]
  list3 = [(i,i=1,0)]
  
  call testsub()
  call testsub([1,2,3])
  call testsub(list1)
  call testsub([integer::])
  call testsub(list2)
  call testsub([(i,i=1,0)])
  call testsub(list3)
  
  call testsub2()
  call testsub2([1,2,3])
  call testsub2(list1)
  call testsub2([integer::])
  call testsub2(list2)
  call testsub2([(i,i=1,0)])
  call testsub2(list3)
end subroutine

subroutine testsub(this)
  implicit none
  
  integer, optional :: this(:)
  
  if (present(this)) then
    call print_line('PRESENT')
  else
    call print_line('not present')
  endif
end subroutine

subroutine testsub2(this)
  implicit none
  
  integer, optional :: this(:)
  
  call testsub(this)
end subroutine
end module

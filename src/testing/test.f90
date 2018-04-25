! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test_keywords
  public :: test_mode
  public :: test
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function test_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function test_mode() result(output)
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
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(String), allocatable :: temp(:)
  type(String) :: two
  
  integer :: i
  
  wd = arguments%value('working_directory')
  
  two = '2'
  
  call print_line(join(left_pad([3,5],'hi')))
  call print_line(join(left_pad([1,2],'2')))
  call print_line(str(maxval([1,2])))
  call print_line(str(2))
  temp = left_pad([1,2],'2')
  temp = left_pad([1,2],two)
  temp = left_pad([1,2],str(2))
  temp = left_pad([1,2],str(maxval([1,2])))
  call print_line(join(left_pad([1,2],str(maxval([1,2])))))
end subroutine
end module

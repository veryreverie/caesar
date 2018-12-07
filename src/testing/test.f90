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
  
  integer :: i
  
  do i=0,4
    call print_line(exp(log_factorial(i)))
  enddo
  
  call print_line('')
  call print_line(exp(log_factorial(4)))
  call print_line(exp(log_factorial(2)))
  call print_line(exp(log_factorial(2)+log_factorial(2)))
  call print_line(exp(log_factorial(4)-log_factorial(2)-log_factorial(2)))
  call print_line(real_multinomial(4,[2,2]))
end subroutine
end module

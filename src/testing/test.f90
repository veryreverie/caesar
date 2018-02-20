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
  use univariate_module
  use monomial_module
  use polynomial_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  type(Univariate) :: u1,u2,u3
  type(Monomial)   :: m1,m2,m3
  type(Polynomial) :: p
  
  integer :: i
  
  wd = arguments%value('working_directory')
  
  u1 = Univariate(1,1)
  u2 = Univariate(2,3)
  u3 = Univariate(3,2)
  call print_line('')
  call print_line('u1^n = '//u1)
  call print_line('u2^n = '//u2)
  call print_line('u3^n = '//u3)
  
  m1 = Monomial(3.14_dp, [u1,u2,u3])
  m2 = Monomial(1.41_dp, [u1,u3])
  m3 = Monomial(2.17_dp, [u2])
  call print_line('')
  call print_line('m1 = '//m1)
  call print_line('m2 = '//m2)
  call print_line('m3 = '//m3)
  
  p = Polynomial([m1,m2,m3])
  call print_line('')
  call print_line('p=')
  call print_line(p)
  
  do i=1,3
    call print_line('')
    call print_line('dp/du'//i)
    call print_line(p%derivative(i))
  enddo
  
end subroutine
end module

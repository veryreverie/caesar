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
  
  logical :: t=.true.
  logical :: f=.false.
  
  logical, allocatable :: a(:)
  logical, allocatable :: m(:)
  
  integer, allocatable :: b(:)
  
  integer :: i
  
  a = [t,f,f,f,f,t]
  call print_line('A:        '//a)
  call print_line('first(a): '//first(a))
  call print_line('last(a) : '//last(a))
  
  a = [f,f,t,f,t,f,f]
  call print_line('A:        '//a)
  call print_line('first(a): '//first(a))
  call print_line('last(a) : '//last(a))
  
  a = [f,f,t,f,f,t,t,t,f,f,t,f,t,t,t,f,f,t,f,t,f,f,t,f,t,f,f]
  m = [f,f,f,t,f,t,t,f,t,f,f,f,t,f,t,f,t,t,f,t,t,f,t,f,t,f,t]
  call print_line('A:          '//a)
  call print_line('first(a)  : '//first(a))
  call print_line('last(a)   : '//last(a))
  call print_line('first(a,m): '//first(a,mask=m))
  call print_line('last(a,m) : '//last(a,mask=m))
  
  a = [f,f,f,f,f,f,f,f,f,f,t,t,t,t,t,t,t,t,t,t,t,t,t,t,t,t,t]
  call print_line('A:          '//a)
  call print_line('first(a)  : '//first(a))
  call print_line('first(a,m): '//first(a,mask=m))
  call print_line('first(a)  : '//first(a,sorted=.true.))
  call print_line('first(a,m): '//first(a,mask=m,sorted=.true.))
  
  a = [t,t,t,t,t,t,t,t,t,t,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f]
  call print_line('A:         '//a)
  call print_line('last(a)  : '//last(a))
  call print_line('last(a,m): '//last(a,mask=m))
  call print_line('last(a)  : '//last(a,sorted=.true.))
  call print_line('last(a,m): '//last(a,mask=m,sorted=.true.))
  
  b = [1,3,4,5,7,9,15,17,19,22,25,28,30]
  call print_line('b       : '//b)
  call print_line('fe(b,25) : '//first_equivalent(b,25,default=0))
  call print_line('fe(b,25) : '//first_equivalent(b,25,sorted=.true.,default=0))
  call print_line('fe(b,24) : '//first_equivalent(b,24,default=0))
  call print_line('fe(b,24) : '//first_equivalent(b,24,sorted=.true.,default=0))
  
contains
  function is_one(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(integer)
      output = input==1
    end select
  end function
end subroutine
end module

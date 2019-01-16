! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  implicit none
  
  private
  
  public :: test
  
  type, abstract :: A
  contains
    procedure(shout_A), public, deferred :: shout
  end type
  
  abstract interface
    subroutine shout_A(this)
      import A
      implicit none
      
      class(A), intent(in) :: this
    end subroutine
  end interface
  
  type, extends(A) :: B
    integer :: b
  contains
    procedure, public :: shout => shout_B
  end type
  
  type, extends(A) :: C
    class(A), allocatable :: c_
    type(String)          :: st
  contains
    procedure, public :: shout => shout_C
  end type
  
  interface C
    module procedure new_C
  end interface
  
contains

subroutine shout_B(this)
  implicit none
  
  class(B), intent(in) :: this
  
  call print_line('B'//this%b)
end subroutine

subroutine shout_C(this)
  implicit none
  
  class(C), intent(in) :: this
  
  call this%c_%shout()
end subroutine

function new_C(input) result(output)
  implicit none
  
  class(A), intent(in) :: input
  type(C)              :: output
  
  select type(input); type is(C)
    output = input
  type is(B)
    allocate(output%c_, source=input)
    output%st = 'c'
  class default
    call err()
  end select
end function

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
  
  class(A), allocatable :: a_
  
  type(C) :: c_
  
  class(A), allocatable :: as_(:)
  class(C), allocatable :: cs_(:)
  
  call print_line('')
  allocate(a_, source=B(1))
  call a_%shout()
  
  call print_line('')
  select type(a_); type is(B)
    call print_line('Type is B')
  end select
  
  call print_line('')
  c_ = C(a_)
  call c_%shout()
  
  call print_line('')
  allocate(as_(4), source=B(1))
  do i=1,size(as_)
    call as_(i)%shout()
  enddo
  
  call print_line('')
  deallocate(as_)
  allocate(as_, source=[B(1),B(2)])
  do i=1,size(as_)
    call as_(i)%shout()
  enddo
  
  call print_line('')
  cs_ = [C(B(1)), C(B(2)), C(B(3))]
end subroutine
end module

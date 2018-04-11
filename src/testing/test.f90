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

subroutine thing(a)
  implicit none
  
  logical, intent(in), optional :: a
  
  if (present(a).lazyand.a) then
    call print_line('True')
  else
    call print_line('Not present / False.')
  endif
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  logical, allocatable :: t
  logical, allocatable :: f
  logical, allocatable :: n(:)
  
  integer :: i
  
  wd = arguments%value('working_directory')
  
  t = .true.
  f = .false.
  
  allocate(n(2))
  n(1) = .true.
  n(2) = .false.
  do i=1,3
    if (i<=size(n) .lazyand. n(i)) then
      call print_line('n('//i//')=.true.')
    else
      call print_line('n('//i//')=.false.')
    endif
  enddo
  
  call thing(.true.)
  call thing(.false.)
  call thing()
end subroutine
end module

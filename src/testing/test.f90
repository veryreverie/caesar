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
  
  type(RandomReal) :: randgen
  
  real(dp) :: random_real
  
  character(8)  :: date
  character(10) :: time
  integer       :: values(8)
  
  integer :: i
  
  wd = arguments%value('working_directory')
  
  do i=1,5
    call random_number(random_real)
    call print_line(random_real)
  enddo
  
  randgen = RandomReal(int(arguments%value('random_seed')))
  
  call print_line('')
  do i=1,5
    call print_line(randgen%random_number())
  enddo
  
  call print_line('')
  do i=1,5
    randgen = RandomReal(randgen%get_seed())
    call print_line(randgen%random_number())
  enddo
  
  do i=1,5
    call print_line(randgen%random_numbers(3))
  enddo
end subroutine
end module

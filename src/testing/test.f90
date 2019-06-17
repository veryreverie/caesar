! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  use anharmonic_module
  
  implicit none
  
  private
  
  public :: startup_test
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_test()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'test'
  mode%description = 'Runs temporary code for testing purposes.'
  mode%keywords = [KeywordData::]
  mode%main_subroutine => test_subroutine
  mode%suppress_from_helptext = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main function.
! ----------------------------------------------------------------------
subroutine test_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  integer, allocatable :: arg
  integer, allocatable :: arg2(:)
  
  call opt(arg)
  call nopt(arg)
  arg = 2
  call opt(arg)
  
  call opt2(arg2)
  arg2 = [integer::]
  call opt2(arg2)
  arg2 = [1]
  call opt2(arg2)
end subroutine

subroutine opt(arg)
  implicit none
  
  integer, intent(in), optional :: arg
  
  if (present(arg)) then
    call print_line('PRESENT')
    call print_line(arg)
  else
    call print_line('NOT PRESENT')
  endif
end subroutine

subroutine nopt(arg)
  implicit none
  
  integer, intent(in) :: arg
  
  call print_line(arg)
end subroutine

subroutine opt2(arg)
  implicit none
  
  integer, intent(in), optional :: arg(:)
  
  if (present(arg)) then
    call print_line('PRESENT')
    call print_line(join(arg))
  else
    call print_line('NOT PRESENT')
  endif
end subroutine
end module

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
  
  integer :: ialloc
  
  mode%mode_name = 'test'
  mode%description = 'Runs temporary code for testing purposes.'
  allocate(mode%keywords(0), stat=ialloc); call err(ialloc)
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
  
  
  logical, allocatable :: a
  
  call print_line(set_default(a,.true.))
  call print_line(set_default(a,.false.))
  a = .true.
  call print_line(set_default(a,.true.))
  call print_line(set_default(a,.false.))
  a = .false.
  call print_line(set_default(a,.true.))
  call print_line(set_default(a,.false.))
end subroutine
end module

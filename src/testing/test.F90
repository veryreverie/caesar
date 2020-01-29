! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  use anharmonic_module
  
  implicit none
  
  private
  
  public :: startup_test
  
  type :: T
  end type
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
  
  integer, allocatable :: b
  type(T), allocatable :: d
  
  type(T) :: a
  
  call print_line(allocated(b))
  call f1(1,b,1)
  
  call print_line(allocated(d))
  call f2(d)
end subroutine

subroutine f1(a,b,c)
  implicit none
  
  integer, intent(in)           :: a
  integer, intent(in), optional :: b
  integer, intent(in)           :: c
  
  call print_line(present(b))
end subroutine

subroutine f2(b)
  implicit none
  
  class(T), intent(in), optional :: b
  
  call print_line(present(b))
end subroutine
end module

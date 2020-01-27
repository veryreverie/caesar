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
  
  call print_line(1 .in. [1,2,3])
  call print_line(1 .in. [2,3])
  call print_line([1,3] .in. [2,3])
  call print_line([4,3] .in. [2,3])
  call print_line([integer::] .in. [2,3])
  call print_line([4,3] .in. [integer::])
  call print_line([integer::] .in. [integer::])
end subroutine
end module

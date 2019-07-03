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
  
  call pulay_solver_example()
end subroutine
end module

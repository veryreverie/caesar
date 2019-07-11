! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  use anharmonic_module
  
  implicit none
  
  private
  
  public :: startup_test
  
  type :: TT
    integer :: x
  end type
  
#define MACRO_TYPE_NAME TT
#include "array_operations.fpp"
contains
#define MACRO_BODY
#include "array_operations.fpp"

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
  
  type(TT), allocatable :: a(:)
  
  a = [TT(1), TT(3)]
  call print_line(a%x)
  
  call append(a, TT(5))
  call print_line(a%x)
  
  call append(a, a)
  call print_line(a%x)
end subroutine
end module

! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module test_module
  use common_module
  
  use anharmonic_module
  
  implicit none
  
  private
  
  public :: startup_test
  
  type :: T1
  contains
    procedure, public :: print => print_T1
  end type
  
  type, extends(T1) :: T2
  contains
    procedure, public :: print => print_T2
  end type
contains

subroutine print_T1(this)
  implicit none
  
  class(T1), intent(in) :: this
  
  call print_line('T1')
end subroutine

subroutine print_T2(this)
  implicit none
  
  class(T2), intent(in) :: this
  
  call print_line('T2')
end subroutine

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
  
  type(T1) :: t_1
  type(T2) :: t_2
  
  call t_1%print()
  call t_2%print()
end subroutine
end module

! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module a_module
  implicit none
  
  private
  
  public :: A
  
  type, abstract :: A
  contains
    generic, public :: temp => &
                     & temp_1
    procedure(temp_1_A), private, deferred :: temp_1
  end type
  
  abstract interface
    subroutine temp_1_A(this)
      import A
      implicit none
      
      class(A), intent(in) :: this
    end subroutine
  end interface
end module

module b_module
  use a_module
  implicit none
  
  private
  
  public :: B
  
  type, extends(A) :: B
  contains
    procedure, private :: temp_1 => temp_1_B
  end type
contains

subroutine temp_1_B(this)
  implicit none
  
  class(B), intent(in) :: this
end subroutine
end module

module test_module
  use common_module
  implicit none
  
  private
  
  public :: test
contains

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
  
  type(String) :: wd
  
  wd = arguments%value('working_directory')
  
end subroutine
end module

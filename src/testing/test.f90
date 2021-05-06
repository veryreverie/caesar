! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module caesar_test_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: test_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module function test_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main module function.
    ! ----------------------------------------------------------------------
    module subroutine test_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

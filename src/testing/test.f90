! ======================================================================
! A test space, for temporary checking of misc. parts of the code.
! ======================================================================
module caesar_test_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  use caesar_potential_data_module
  use caesar_anharmonic_module
  implicit none
  
  private
  
  public :: startup_test
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_test() 
    end subroutine
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

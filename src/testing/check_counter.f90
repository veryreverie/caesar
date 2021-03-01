! ======================================================================
! Checks SharedCounter to see whether a gfortran bug is present or not.
! ======================================================================
module caesar_check_counter_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_check_counter
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_check_counter() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main module function.
    ! ----------------------------------------------------------------------
    module subroutine check_counter_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
  
  interface
    module subroutine check_counter_subroutine_2(a) 
      type(SharedCounter), intent(inout) :: a
    end subroutine
  end interface
  
  interface
    module function colour_check(input,expected) result(output) 
      logical, intent(in) :: input
      logical, intent(in) :: expected
      type(String)        :: output
    end function
  end interface
end module

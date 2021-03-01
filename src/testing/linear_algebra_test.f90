! ======================================================================
! Runs tests on linear_algebra.
! ======================================================================
module caesar_linear_algebra_test_module
  use caesar_common_module
  implicit none
  
  interface
    module subroutine linear_algebra_test(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module

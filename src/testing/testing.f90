! ======================================================================
! Provides tests.
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module caesar_testing_module
  use caesar_common_module
  
  use caesar_check_counter_module
  use caesar_update_basis_functions_module
  use caesar_test_module
  implicit none
  
  interface
    module function testing_modes() result(output)
      type(ProgramMode), allocatable :: output(:)
    end function
  end interface
end module

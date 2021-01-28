! ======================================================================
! Provides tests.
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module caesar_testing_module
  use caesar_check_counter_module
  use caesar_update_basis_functions_module
  use caesar_test_module
  implicit none
contains
subroutine startup_testing()
  implicit none
  
  call startup_check_counter()
  call startup_update_basis_functions()
  call startup_test()
end subroutine
end module

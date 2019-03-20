! ======================================================================
! Provides tests.
! ======================================================================
! This module is simply an interface to various subsidiary modules.
module testing_module
  use check_counter_module
  use update_basis_functions_module
  use test_module
  implicit none
contains
subroutine startup_testing()
  implicit none
  
  call startup_check_counter()
  call startup_update_basis_functions()
  call startup_test()
end subroutine
end module

! ======================================================================
! A dummy module for when QUIP is not being linked against.
! ======================================================================
module quip_wrapper_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

subroutine test_quip()
  implicit none
  
  call err()
end subroutine
end module

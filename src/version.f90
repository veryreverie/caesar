! ======================================================================
! Prints the version number.
! ======================================================================
module version_module
  use common_module
  implicit none
contains
subroutine print_version()
  implicit none
  
  call print_line('Caesar version: 0.0.190528A (Gaius)')
end subroutine
end module

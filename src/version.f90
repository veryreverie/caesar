! ======================================================================
! Prints the version number.
! ======================================================================
module caesar_version_module
  use caesar_common_module
  implicit none
contains
subroutine print_version()
  implicit none
  
  call print_line('Caesar version: 0.1.210319C (Julius)')
end subroutine
end module

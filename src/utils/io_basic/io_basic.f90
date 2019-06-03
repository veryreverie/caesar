! ======================================================================
! Basic I/O functionality, including:
!    - Variable-length strings.
!    - Formatted and error-checked writing to the terminal.
!    - System calls.
! ======================================================================
! This module is simply an interface for the various I/O submodules.
module io_basic_module
  use error_module
  use string_module
  use print_settings_module
  use print_module
  use intrinsics_module
  use io_utils_module
  use token_module
  implicit none
contains
subroutine startup_io_basic()
  implicit none
  
  call startup_io_utils()
end subroutine
end module

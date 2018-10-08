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
  implicit none
end module

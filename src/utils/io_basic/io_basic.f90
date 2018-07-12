! ======================================================================
! Basic I/O functionality, including:
!    - Variable-length strings.
!    - Formatted and error-checked writing to the terminal.
!    - System calls.
! ======================================================================
! This module is simply an interface for the various I/O submodules.
module io_basic_module
  use error_submodule
  use string_submodule
  use print_settings_submodule
  use print_submodule
  use intrinsics_submodule
  use io_utils_submodule
  implicit none
end module

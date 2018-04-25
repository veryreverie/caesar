! ======================================================================
! I/O functionality, including:
!    - Variable-length strings.
!    - Formatted and error-checked writing to the terminal.
!    - File handling.
!    - Base types from which all extended types can be printed out.
!    - System calls.
! ======================================================================
! This module is simply an interface for the various I/O submodules.
module io_module
  use error_submodule
  use string_submodule
  use print_submodule
  use intrinsics_submodule
  use stringable_submodule
  use printable_submodule
  use io_submodule
  use ifile_submodule
  use ofile_submodule
  implicit none
end module

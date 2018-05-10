! ======================================================================
! More advanced I/O functionality, including:
!    - File handling.
!    - Base types from which extended types can be written and read.
! ======================================================================
! This module is simply an interface for the various I/O submodules.
module io_module
  use io_basic_module
  
  use string_writeable_submodule
  use strings_writeable_submodule
  use string_array_submodule
  use ifile_submodule
  use ofile_submodule
  use string_readable_submodule
  use strings_readable_submodule
  use stringable_submodule
  use stringsable_submodule
  implicit none
end module

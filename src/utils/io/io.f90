! ======================================================================
! More advanced I/O functionality, including:
!    - File handling.
!    - Base types from which extended types can be written and read.
! ======================================================================
! This module is simply an interface for the various I/O submodules.
module io_module
  use io_basic_module
  
  use string_writeable_module
  use strings_writeable_module
  use string_array_module
  use ifile_module
  use ofile_module
  use string_readable_module
  use strings_readable_module
  use stringable_module
  use stringsable_module
  implicit none
end module

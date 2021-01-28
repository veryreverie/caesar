! ======================================================================
! More advanced I/O functionality, including:
!    - File handling.
!    - Base types from which extended types can be written and read.
! ======================================================================
! This module is simply an interface for the various I/O modules.
module caesar_io_module
  use caesar_io_basic_module
  
  use caesar_string_writeable_module
  use caesar_strings_writeable_module
  use caesar_string_array_module
  use caesar_ifile_module
  use caesar_ofile_module
  use caesar_string_readable_module
  use caesar_strings_readable_module
  use caesar_stringable_module
  use caesar_stringsable_module
  implicit none
contains
subroutine startup_io()
  implicit none
  
  call startup_io_basic()
end subroutine
end module

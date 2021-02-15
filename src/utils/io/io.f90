!> This module is an interface for the various I/O modules.
module caesar_io_module
  use caesar_error_module
  use caesar_string_base_module
  use caesar_string_module
  use caesar_print_settings_module
  use caesar_print_module
  use caesar_intrinsics_module
  use caesar_io_utils_module
  use caesar_token_module
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
  
  call startup_io_utils()
end subroutine
end module

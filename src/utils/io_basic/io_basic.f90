!> This module is an interface for the various I/O modules.
module caesar_io_basic_module
  use caesar_error_module
  use caesar_string_base_module
  use caesar_string_module
  use caesar_print_settings_module
  use caesar_print_module
  use caesar_intrinsics_module
  use caesar_io_utils_module
  use caesar_token_module
  implicit none
contains
subroutine startup_io_basic()
  implicit none
  
  call startup_io_utils()
end subroutine
end module

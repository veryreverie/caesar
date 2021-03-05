submodule (caesar_current_working_directory_module) caesar_current_working_directory_submodule
  use caesar_print_module
  use caesar_error_module
  use caesar_c_string_module
  
  logical      :: INITIALISED = .false.
  type(String) :: STORED_CURRENT_WORKING_DIRECTORY
contains

!> Helper function for [[current_working_directory]].
!> Gets the current working directory `.`
!>    by calling [[get_current_working_directory_c]].
function get_current_working_directory() result(output)
  type(String) :: output
  
  integer, parameter :: result_size = 1024
  
  character(result_size) :: current_dir
  logical                :: success
  
  success = get_current_working_directory_c(result_size, current_dir)
  
  if (.not. success) then
    call print_line(ERROR//': getcwd failed.')
    call err()
  endif
  
  output = parse_c_string(current_dir)
end function

module procedure current_working_directory
  if (.not. INITIALISED) then
    STORED_CURRENT_WORKING_DIRECTORY = get_current_working_directory()
    INITIALISED = .true.
  endif
  
  output = STORED_CURRENT_WORKING_DIRECTORY
end procedure

module procedure set_current_working_directory_String
  STORED_CURRENT_WORKING_DIRECTORY = working_directory
  INITIALISED = .true.
end procedure

module procedure set_current_working_directory_character
  STORED_CURRENT_WORKING_DIRECTORY = working_directory
  INITIALISED = .true.
end procedure
end submodule

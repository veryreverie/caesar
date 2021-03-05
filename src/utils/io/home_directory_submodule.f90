submodule (caesar_home_directory_module) caesar_home_directory_submodule
  use caesar_print_module
  use caesar_error_module
  use caesar_c_string_module
  
  logical      :: INITIALISED = .false.
  type(String) :: STORED_HOME_DIRECTORY
contains

!> Helper function for [[home_directory]].
!> Gets the home `~` directory by calling [[get_home_directory_c]].
function get_home_directory() result(output)
  type(String) :: output
  
  integer, parameter  :: result_size = 1024
  
  character(result_size) :: home_dir
  logical                :: success
  
  success = get_home_directory_c(home_dir)
  
  if (.not. success) then
    call print_line(ERROR//': getenv("HOME") failed.')
    call err()
  endif
  
  output = parse_c_string(home_dir)
end function

module procedure home_directory
  if (.not. INITIALISED) then
    STORED_HOME_DIRECTORY = get_home_directory()
    INITIALISED = .true.
  endif
  
  output = STORED_HOME_DIRECTORY
end procedure
end submodule

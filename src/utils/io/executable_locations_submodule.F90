submodule (caesar_executable_locations_module) caesar_executable_locations_submodule
  use caesar_print_module
  use caesar_error_module
  use caesar_c_string_module
  
  logical      :: INITIALISED = .false.
  type(String) :: STORED_EXECUTABLE_LOCATION
  type(String) :: STORED_PYTHON_SCRIPTS_LOCATION
contains

subroutine initialise()
  if (.not. INITIALISED) then
    ! EXECUTABLE_LOCATION and PYTHON_SCRIPTS_LOCATION are
    !    preprocessor variables.
    STORED_EXECUTABLE_LOCATION = EXECUTABLE_LOCATION
    STORED_PYTHON_SCRIPTS_LOCATION = PYTHON_SCRIPTS_LOCATION
    INITIALISED = .true.
  endif
end subroutine

module procedure executable_location
  call initialise()
  
  output = STORED_EXECUTABLE_LOCATION
end procedure

module procedure python_scripts_location
  call initialise()
  
  output = STORED_PYTHON_SCRIPTS_LOCATION
end procedure
end submodule

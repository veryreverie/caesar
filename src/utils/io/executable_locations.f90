!> Provides the [[executable_location]] and [[pyton_scripts_location]]
!>    functions.
module caesar_executable_locations_module
  use caesar_string_base_module
  use caesar_string_module
  implicit none
  
  private
  
  public :: executable_location
  public :: python_scripts_location
  
  interface
    !> Returns the path to the `caesar` executable location.
    module function executable_location() result(output) 
      type(String) :: output
    end function
    
    !> Returns the path to the directory containing Caesar's python scripts.
    module function python_scripts_location() result(output)
      type(String) :: output
    end function
  end interface
end module

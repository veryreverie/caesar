!> Provides the [[home_directory]] function.
module caesar_home_directory_module
  use caesar_string_base_module
  use caesar_string_module
  implicit none
  
  private
  
  public :: home_directory
  
  interface
    !> Gets the home `~` directory by calling `getenv('HOME')`.
    function get_home_directory_c(home) bind(c) result(success)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), intent(out) :: home(*)
      logical(kind=c_bool)                :: success
    end function
  
    !> Returns the path to the home `~` directory.
    module function home_directory() result(output)
      type(String) :: output
    end function
  end interface
end module

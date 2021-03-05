!> Provides the [[parse_c_string]] function.
module caesar_c_string_module
  use caesar_string_base_module
  use caesar_string_module
  implicit none
  
  private
  public :: parse_c_string

  interface parse_c_string
    !> Converts a `NUL`-terminated C string to a Fortran [[String(type)]].
    module function parse_c_string(input) result(output) 
      character(*), intent(in) :: input
      type(String)             :: output
    end function
  end interface
end module

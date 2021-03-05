!> Provides the [[current_working_directory]] and
!>    [[set_current_working_directory]] functions.
module caesar_current_working_directory_module
  use caesar_string_base_module
  use caesar_string_module
  implicit none
  
  private
  
  public :: current_working_directory
  public :: set_current_working_directory
  
  interface
    !> Gets the current working directory `.` by calling `getcwd`.
    function get_current_working_directory_c(result_size, cwd) bind(c) &
       & result(success)
      use, intrinsic :: iso_c_binding
      integer(kind=c_int),    intent(in)  :: result_size
      character(kind=c_char), intent(out) :: cwd(*)
      logical(kind=c_bool)                :: success
    end function
    
    !> Returns the path to the current working directory `.`.
    module function current_working_directory() result(output) 
      type(String) :: output
    end function
  end interface
  
  interface set_current_working_directory
    !> Sets the current working directory to `working_directory`.
    module subroutine set_current_working_directory_String(working_directory) 
      type(String), intent(in) :: working_directory
    end subroutine
    
    !> Sets the current working directory to `working_directory`.
    module subroutine set_current_working_directory_character( &
       & working_directory) 
      character(*), intent(in) :: working_directory
    end subroutine
  end interface
end module

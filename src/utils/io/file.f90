!> Provides file handling methods which are common to the [[IFile(type)]] and
!>    [[OFile(type)]] classes.
!> Should not be imported by any modules other than ifile and ofile.
module caesar_file_module
  use caesar_foundations_module
  implicit none
  
  private
  
  public :: open_read_file
  public :: open_write_file
  
  interface open_file
    !> Open a file, and return the unit it is opened in.
    !> See the `open` intrinsic for details on `status`, `action` and `access`.
    module function open_file(filename,status,action,access) result(unit_num) 
      character(*), intent(in) :: filename
      character(*), intent(in) :: status
      character(*), intent(in) :: action
      character(*), intent(in) :: access
      integer                  :: unit_num
    end function
  end interface
  
  interface open_read_file
    !> Open a file for reading, and return the unit it is opened in.
    module function open_read_file(filename) result(unit_num) 
      character(*), intent(in) :: filename
      integer                  :: unit_num
    end function
  end interface
  
  interface open_write_file
    !> Open a file for writing, and return the unit it is opened in.
    module function open_write_file(filename) result(unit_num) 
      character(*), intent(in) :: filename
      integer                  :: unit_num
    end function
  end interface
end module

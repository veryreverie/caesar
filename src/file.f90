module file_io
  use string_module
  implicit none
  
  private
  
  public :: open_read_file   ! open a file for reading
  public :: open_write_file  ! open a file for writing
  public :: open_append_file ! open a file for appending to
  public :: file_exists      ! checks if a file exists
  public :: count_lines      ! counts the number of lines in a file
  
  interface open_read_file
    module procedure open_read_file_character
    module procedure open_read_file_string
  end interface
  
  interface open_write_file
    module procedure open_write_file_character
    module procedure open_write_file_string
  end interface
  
  interface open_append_file
    module procedure open_append_file_character
    module procedure open_append_file_string
  end interface
  
  interface file_exists
    module procedure file_exists_character
    module procedure file_exists_string
  end interface

contains

! open a file with a specified mode, and return the unit it is opened in
function open_file(filename,status,action,access) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  character(*), intent(in) :: status
  character(*), intent(in) :: action
  character(*), intent(in) :: access
  integer                  :: unit_num
  
  integer :: iostat
  logical :: opened
  
  unit_num = 0
  opened=.true.
  
  do while (opened) ! loop until unused unit is found
    unit_num = unit_num + 1
    if (unit_num==5 .or. unit_num==6) cycle ! ignore 0,5,6 std units
    if (unit_num>=100 .and. unit_num<=102) cycle ! ignore 100,101,102 std units
    
    inquire(unit=unit_num,opened=opened,iostat=iostat) ! check unit
    if (iostat/=0) opened=.true. ! ignore unit if error
  enddo
  
  open(unit=unit_num,file=filename,status=status,action=action,&
    &access=access,iostat=iostat)
  if (iostat /= 0) then
    write(*,*) 'Problem opening '//trim(filename)//' file.'
    stop
  endif
end function

! open a file for reading, and return the unit it is opened in
function open_read_file_character(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'old', 'read', 'sequential')
end function

function open_read_file_string(filename) result(unit_num)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_read_file(char(filename))
end function

! open a file for writing, and return the unit it is opened in
function open_write_file_character(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'sequential')
end function

function open_write_file_string(filename) result(unit_num)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_write_file(char(filename))
end function

! open a file for appending, and return the unit it is opened in
function open_append_file_character(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'append')
end function

function open_append_file_string(filename) result(unit_num)
  implicit none
  
  type(String), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_append_file(char(filename))
end function

! ----------------------------------------------------------------------
! Checks if a file exists
! ----------------------------------------------------------------------
function file_exists_character(filename) result(output)
  implicit none
  
  character(*), intent(in) :: filename
  logical                  :: output
  
  inquire(file=filename, exist=output)
end function

function file_exists_string(filename) result(output)
  implicit none
  
  type(String), intent(in) :: filename
  logical                  :: output
  
  output = file_exists(char(filename))
end function

! ----------------------------------------------------------------------
! Gets the number of lines remaining in a file
! Rewinds the file back to beginning
! ----------------------------------------------------------------------
function count_lines(file_unit) result(output)
  implicit none
  
  integer, intent(in) :: file_unit
  integer             :: output
  
  integer       :: eof_reached
  character(80) :: line
  
  output = 0
  do while (eof_reached==0)
    read(file_unit, *, iostat=eof_reached) line
    if (eof_reached==0) then
      output = output+1
    endif
  enddo
  rewind(file_unit)
end function

end module

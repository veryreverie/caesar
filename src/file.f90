module file_io
  implicit none
  
  private :: open_file        ! helper function
  public  :: open_read_file   ! open a file for reading
  public  :: open_write_file  ! open a file for writing
  public  :: open_append_file ! open a file for appending to
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
function open_read_file(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'old', 'read', 'sequential')
end function

! open a file for writing, and return the unit it is opened in
function open_write_file(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'sequential')
end function

! open a file for appending, and return the unit it is opened in
function open_append_file(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'append')
end function
end module

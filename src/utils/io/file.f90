! ======================================================================
! File handling common to both input and output files.
! Should not be imported by any modules other than ifile and ofile.
! ======================================================================
module file_submodule
  use precision_module
  use io_basic_module
  implicit none
  
  private
  
  public :: open_read_file
  public :: open_write_file
contains

! ----------------------------------------------------------------------
! open a file with a specified mode, and return the unit it is opened in.
! ----------------------------------------------------------------------
function open_file(filename,status,action,access) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  character(*), intent(in) :: status
  character(*), intent(in) :: action
  character(*), intent(in) :: access
  integer                  :: unit_num
  
  type(String) :: formatted_filename
  
  integer :: iostat
  logical :: opened
  
  logical :: unit_found
  
  unit_num = 0
  unit_found = .false.
  
  ! Loop over file units until an unused unit is found.
  do while (.not. unit_found)
    unit_num = unit_num + 1
    
    ! Ignore standard units 0, 5 and 6.
    if (unit_num==5 .or. unit_num==6) then
      cycle
    endif
    
    ! Ignore standard units 100, 101 and 102.
    if (unit_num>=100 .and. unit_num<=102) then
      cycle
    endif
    
    ! Check file unit.
    inquire(unit=unit_num,opened=opened,iostat=iostat)
    if (iostat==0 .and. .not. opened) then
      unit_found = .true.
    endif
  enddo
  
  ! Format filename.
  formatted_filename = format_path(filename)
  
  open( unit   = unit_num,                 &
      & file   = char(formatted_filename), &
      & status = status,                   &
      & action = action,                   &
      & access = access,                   &
      & iostat = iostat                    )
  if (iostat /= 0) then
    call print_line('Error opening '//formatted_filename//' file.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! open a file for reading, and return the unit it is opened in.
! ----------------------------------------------------------------------
function open_read_file(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'old', 'read', 'sequential')
end function

! ----------------------------------------------------------------------
! open a file for writing, and return the unit it is opened in.
! ----------------------------------------------------------------------
function open_write_file(filename) result(unit_num)
  implicit none
  
  character(*), intent(in) :: filename
  integer                  :: unit_num
  
  unit_num = open_file(filename, 'unknown', 'write', 'sequential')
end function
end module

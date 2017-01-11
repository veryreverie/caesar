module process
  implicit none
  
  ! holds the status, stderr and stdout of a system process
  type ProcessResult
    integer        :: status
    character(100) :: stdout
    character(100) :: stderr
  end type
  
contains

! calls system(input), returning a ProcessResult
function system_process(input) result(output)
  use file_io, only : open_read_file
  implicit none
  
  character(len=*), intent(in) :: input
  type(ProcessResult)          :: output
  
  character(100) :: temp_char
  integer        :: file_unit
  integer        :: iostat
  
  ! run process
  output%status = system(input//' >system_process_o 2>system_process_e')
  
  ! read stdout
  file_unit = open_read_file('system_process_o')
  read(file_unit,*,iostat=iostat) temp_char
  close(file_unit)
  call system('rm system_process_o')
  if (iostat==0) then
    output%stdout=temp_char
  else
    output%stdout=''
  endif
  
  ! read sterr
  file_unit = open_read_file('system_process_e')
  read(file_unit,*,iostat=iostat) temp_char
  close(file_unit)
  call system('rm system_process_e')
  if (iostat==0) then
    output%stderr=temp_char
  else
    output%stderr=''
  endif
end function
end module

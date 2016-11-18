module process
  implicit none
  
  ! holds the status, stderr and stdout of a system process
  type ProcessResult
    integer       :: status
    character(32) :: stdout
    character(32) :: stderr
  end type
  
contains

! calls system(input), returning a ProcessResult
function system_process(input) result(output)
  implicit none
  
  character(len=*), intent(in) :: input
  type(ProcessResult)          :: output
  
  character(32)                :: temp_char
  integer                      :: io_stat
  
  ! run process
  output%status = system(input//' >system_process_o 2>system_process_e')
  
  ! read stdout
  open(unit=2, file='system_process_o')
  read(2,*,iostat=io_stat) temp_char
  if (io_stat==0) then
    output%stdout = temp_char
  else
    output%stdout = ''
  endif
  close(2)
  call system('rm system_process_o')
  
  ! read sterr
  open(unit=2, file='system_process_e')
  read(2,*,iostat=io_stat) temp_char
  if (io_stat==0) then
    output%stderr = temp_char
  else
    output%stderr = ''
  endif
  close(2)
  call system('rm system_process_e')
end function
end module

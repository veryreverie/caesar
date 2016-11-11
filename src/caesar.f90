module utils
contains

  ! returns an array containing the command line arguments.
  function command_line_args() result(o_args)
    implicit none

    integer                        :: i
    integer                        :: arg_count
    character(len=32), allocatable :: o_args(:)
    
    ! read in the number of arguments
    arg_count = iargc()
    ! allocate the return array
    allocate (o_args(arg_count))
    ! read the arguments into the array
    do i=1,arg_count
      call getarg(i, o_args(i))
    enddo
    
    return
  end function

end module

program caesar
  use utils, only : command_line_args
  implicit none
  
  integer           :: i
  character(len=32), allocatable :: args(:)

  ! read in command line arguments
  args = command_line_args()
  
  if (size(args) == 0) then
    write(*,*) "No arguments given. For help, call caesar -h"
    stop
  elseif (args(1) == "-h" .or. args(1) == "--help") then
    write(*,*) "[Placeholder help text]"
    stop
  else
    write(*,*) "Unrecognised argument : "//args(1)
    stop
  endif
  
  deallocate (args)

end program

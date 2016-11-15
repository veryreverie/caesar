module utils
contains
  
  ! returns an array containing the command line arguments.
  function command_line_args() result(o_args)
    implicit none

    integer                        :: i         ! loop index
    integer                        :: arg_count ! no. of command line args
    character(len=32), allocatable :: o_args(:) ! return value

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

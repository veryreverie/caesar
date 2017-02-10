! ======================================================================
! various utilities
! ======================================================================
module utils
  use constants,      only : dp
  implicit none
  
  interface lower_case
    module procedure lower_case_character
  end interface
  
contains

! ----------------------------------------------------------------------
! reports an error and stops
! ----------------------------------------------------------------------
subroutine errstop(subroutine_name, message)
  implicit none
  
  character(*), intent(in) :: subroutine_name ! where errstop is called
  character(*), intent(in) :: message         ! error message
  
  write(*,*) ''
  write(*,*) 'error in subroutine '//trim(adjustl(subroutine_name))//'.'
  write(*,*) ''
  call wordwrap(trim(adjustl(message)))
  write(*,*) ''
  stop
end subroutine

! ----------------------------------------------------------------------
! reports an allocation error and stops
! ----------------------------------------------------------------------
subroutine erralloc(arg)
  implicit none
  
  character(*), intent(in) :: arg
  
  write(*,*) ''
  write(*,*) 'Problem allocating '//trim(adjustl(arg))//' array.'
  write(*,*)
end subroutine

! ----------------------------------------------------------------------
! prints out the contents of the character string 'text',
! ensuring that line breaks only occur at space characters. The output
! is written to unit unit_in if this parameter is supplies; otherwise the
! output is written to unit o. The maximum length of each line is given
! by linelength_in if this is supplied; otherwise the default line length
! is 79 characters.
! ----------------------------------------------------------------------
subroutine wordwrap(text, unit_in, linelength_in)
  implicit none
  
  integer, intent(in), optional :: unit_in
  integer, intent(in), optional :: linelength_in
  character(*), intent(in)      :: text
  character(260)                :: temp
  integer                       :: i
  integer                       :: unit
  integer                       :: lentext
  integer                       :: startpoint
  integer                       :: stoppoint
  integer                       :: lastpos
  integer                       :: linelength
  
  ! check if unit_in is supplied
  if (present(unit_in)) then
    unit = unit_in
  else
    unit=6
  endif
  
  lentext = len(trim(text))
  if (lentext<1) then ! text is empty
    write(unit,*) ""
    return
  endif
  
  ! check if linelength_in is supplied
  if (present(linelength_in)) then
    if (linelength_in>=2) then
      linelength=linelength_in
    else
      linelength=2
    endif
  else
    linelength=79
  endif
  
  startpoint=1
  do i=1, huge(1) ! loop over lines
    stoppoint = startpoint+linelength-1
    if (stoppoint<=lentext) then
      lastpos = index(trim(text(startpoint:stoppoint)),' ',.true.)
      if (lastpos>0) stoppoint = startpoint+lastpos-1
    else
      stoppoint = lentext
    endif
    
    if (i==1) then
      ! allow the user to indent the first line, if they wish
      temp = text(startpoint:stoppoint) ! or pathscale.f90 fails to compile
      write(unit,*) trim(temp)
    else
      temp = text(startpoint:stoppoint) ! or pathscale.f90 fails to compile
      write(unit,*) trim(adjustl(temp))
    endif
    
    if (stoppoint==lentext) then
      exit
    else
      startpoint = stoppoint+1
    endif ! finished text?
  enddo
end subroutine

! ----------------------------------------------------------------------
! converts integers to left justified strings that can be printed in the
! middle of a sentence without introducing large amounts of white space.
! ----------------------------------------------------------------------
function i2s(n) result(output)
  integer, intent(in) :: n
  character(12)       :: output
  write(output,"(I0)") n
end function

! ----------------------------------------------------------------------
! This is how i2s used to be written.
! ----------------------------------------------------------------------
!character(12) function i2s(n)
!  implicit none
!  
!  integer, intent(in) :: n                 ! input integer
!  integer             :: i                 ! characters left to process
!  integer             :: j                 ! loop counter
!  integer, parameter  :: ichar0=ichar('0') ! the character '0'
!  
!  i2s = ''
!  i=abs(n)
!  do j=len(i2s),1,-1
!    i2s(j:j)=achar(ichar0+mod(i,10))
!    i=i/10
!    if (i==0) exit
!  enddo
!  if (n<0) then
!    i2s='-'//adjustl(i2s)
!  else
!    i2s=adjustl(i2s)
!  endif
!end function

! ----------------------------------------------------------------------
! Returns an array containing the command line arguments
! ----------------------------------------------------------------------
function command_line_args() result(args)
  use string_module
  implicit none

  integer                   :: i         ! loop index
  integer                   :: arg_count ! no. of command line args
  type(String), allocatable :: args(:)   ! return value
  
  character(1000) :: temp_char

  ! read in the number of arguments
  arg_count = iargc()
  ! allocate the return array
  allocate (args(arg_count))
  ! read the arguments into the array
  do i=1,arg_count
    call getarg(i, temp_char)
    args(i) = trim(temp_char)
  enddo

  return
end function

! ----------------------------------------------------------------------
! Converts a string to lower case
! ----------------------------------------------------------------------
pure function lower_case_character(input) result(output)
  implicit none
  
  character(*), intent(in) :: input
  character(len(input))    :: output
  
  character(*), parameter :: lower_chars = "abcdefghijklmnopqrstuvwxyz"
  character(*), parameter :: upper_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  
  integer :: i,j
  
  output = input
  
  do i=1,len(output)
    j = index(upper_chars, output(i:i))
    if (j/=0) then
      output(i:i) = lower_chars(j:j)
    endif
  enddo
end function

! ----------------------------------------------------------------------
! Makes all values lie in [-0.5-tol,0.5-tol], by adding or subtracting integers
! ----------------------------------------------------------------------
elemental function reduce_interval(input,tol) result(output)
  use constants, only : dp
  implicit none
  
  real(dp), intent(in) :: input
  real(dp), intent(in) :: tol
  real(dp)             :: output
  
  output = modulo(input+0.5_dp+tol,1.0_dp)-0.5_dp-tol
end function

! ----------------------------------------------------------------------
! As above, but for vec*grid -> [-0.5,0.5)*grid
! ----------------------------------------------------------------------
pure function reduce_to_ibz(input,grid) result(output)
  implicit none
  
  integer, intent(in) :: input(3)
  integer, intent(in) :: grid(3)
  integer             :: output(3)
  
  output = modulo(input+grid/2,grid)-grid/2
end function

! ----------------------------------------------------------------------
! Returns true if input is less than tol from an integer
! ----------------------------------------------------------------------
elemental function is_int(input,tol) result(output)
  use constants, only : dp
  implicit none
  
  real(dp), intent(in) :: input
  real(dp), intent(in) :: tol
  logical              :: output
  
  output = abs(nint(input)-input) < tol
end function
end module

! ======================================================================
! various utilities
! ======================================================================
module utils
  use constants,      only : dp
  implicit none
  
  interface mkdir
    module procedure mkdir_character
    module procedure mkdir_String
  end interface
  
contains

! ----------------------------------------------------------------------
! Reports an error and stops.
! ----------------------------------------------------------------------
subroutine errstop(subroutine_name, message)
  use string_module
  use err_module
  implicit none
  
  character(*), intent(in) :: subroutine_name ! where errstop is called
  character(*), intent(in) :: message         ! error message
  
  call print_line('')
  call print_line('error in subroutine '//trim(adjustl(subroutine_name))//'.')
  call print_line('')
  call print_line(trim(adjustl(message)))
  call print_line('')
  call err()
end subroutine

! ----------------------------------------------------------------------
! reports an allocation error and stops
! ----------------------------------------------------------------------
subroutine erralloc(arg)
  use string_module
  implicit none
  
  character(*), intent(in) :: arg
  
  call print_line('')
  call print_line('Problem allocating '//trim(adjustl(arg))//' array.')
  call print_line('')
end subroutine

! ----------------------------------------------------------------------
! converts integers to left justified strings that can be printed in the
! middle of a sentence without introducing large amounts of white space.
! ----------------------------------------------------------------------
function i2s(n) result(output)
  use string_module
  implicit none
  
  integer, intent(in) :: n
  character(12)       :: output
  
  output = char(str(n))
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

! ----------------------------------------------------------------------
! Make a directory, if it doesn't already exist.
! ----------------------------------------------------------------------
subroutine mkdir_character(dirname)
  use string_module
  implicit none
  
  character(*), intent(in) :: dirname
  
  call system('if [ ! -e '//dirname//' ]; then mkdir '//dirname//'; fi')
end subroutine

subroutine mkdir_String(dirname)
  use string_module
  implicit none
  
  type(String), intent(in) :: dirname
  
  call mkdir(char(dirname))
end subroutine
end module

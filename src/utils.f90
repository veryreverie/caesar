! ======================================================================
! various utilities
! ======================================================================
module utils
  use constants,      only : dp
  implicit none
  
  interface format_path
    module procedure format_path_character
    module procedure format_path_String
  end interface
  
contains

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
  arg_count = command_argument_count()
  ! allocate the return array
  allocate (args(arg_count))
  ! read the arguments into the array
  do i=1,arg_count
    call get_command_argument(i, temp_char)
    args(i) = trim(temp_char)
  enddo

  return
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
! Make a directory, if it doesn't already exist.
! ----------------------------------------------------------------------
subroutine mkdir(dirname)
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: dirname
  
  call system_call('if [ ! -e '//dirname//' ]; then mkdir '//dirname//'; fi')
end subroutine

! ----------------------------------------------------------------------
! Takes a directory name, and converts it into an absolute path
!    in standard format (without a trailing '/').
! ----------------------------------------------------------------------
function format_path_character(path,cwd) result(output)
  use string_module
  use file_module
  implicit none
  
  character(*), intent(in) :: path
  type(String), intent(in) :: cwd       ! Current working directory.
  type(String)             :: output
  
  integer :: last
  
  last = len(path)
  
  if (last==0) then
    call print_line('Error: no path provided.')
    call err()
  endif
  
  ! Trim trailing '/', if present.
  if (path(last:)=='/') then
    last = last - 1
  endif
  
  if (path(:1)=='/' .or. path(:1)=='~') then
    ! Path is absolute.
    output = path(:last)
  else
    ! Path is relative. Prepend current working directory.
    output = cwd//'/'//path(:last)
  endif
end function

function format_path_String(path,cwd) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: path
  type(String), intent(in) :: cwd    ! Current working directory.
  type(String)             :: output
  
  output = format_path(char(path),cwd)
end function

! ----------------------------------------------------------------------
! Converts a file seedname into the appropriate dft input or output filename.
! ----------------------------------------------------------------------
function make_dft_input_filename(dft_code,seedname) result(output)
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (dft_code == 'castep') then
    output = seedname//'.cell'
  elseif (dft_code == 'vasp') then
    output = 'POSCAR'
  elseif (dft_code == 'qe') then
    output = seedname//'.in'
  else
    call print_line('Unrecognised dft code: '//dft_code)
    call err()
  endif
end function

function make_dft_output_filename(dft_code,seedname) result(output)
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: dft_code
  type(String), intent(in) :: seedname
  type(String)             :: output
  
  if (dft_code == 'castep') then
    output = seedname//'.castep'
  elseif (dft_code == 'vasp') then
    output = 'OUTCAR'
  elseif (dft_code == 'qe') then
    output = seedname//'.out'
  else
    call print_line('Unrecognised dft code: '//dft_code)
    call err()
  endif
end function

! ----------------------------------------------------------------------
! Replicates norm2() from f2008 standard.
! ----------------------------------------------------------------------
function l2_norm(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  real(dp)             :: output
  
  output = dsqrt(dot_product(input,input))
end function
end module

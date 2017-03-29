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
  
  interface format_directory
    module procedure format_directory_character
    module procedure format_directory_String
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

! ----------------------------------------------------------------------
! Takes a directory name, and converts it into an absolute path
!    in standard format (without a trailing '/').
! ----------------------------------------------------------------------
function format_directory_character(directory,cwd) result(output)
  use string_module
  use err_module
  implicit none
  
  character(*), intent(in) :: directory
  type(String), intent(in) :: cwd       ! Current working directory.
  type(String)             :: output
  
  integer :: last
  
  last = len(directory)
  
  if (last==0) then
    call print_line('Error: no directory provided.')
    call err()
  endif
  
  ! Trim trailing '/', if present.
  if (directory(last:)=='/') then
    last = last - 1
  endif
  
  if (directory(:1)=='.') then
    ! Directory is relative. Prepend current working directory.
    output = cwd//'/'//directory(:last)
  elseif (directory(:1)=='/') then
    ! Directory is absolute.
    output = directory(:last)
  else
    ! Directory is not valid.
    call print_line('Error: bad directory provided:')
    call print_line(directory)
    call err()
  endif
end function

function format_directory_String(directory,cwd) result(output)
  use string_module
  implicit none
  
  type(String), intent(in) :: directory
  type(String), intent(in) :: cwd       ! Current working directory.
  type(String)             :: output
  
  output = format_directory(char(directory),cwd)
end function

! ----------------------------------------------------------------------
! Converts a file seedname into the appropriate dft input or output filename.
! ----------------------------------------------------------------------
function make_dft_input_filename(dft_code,seedname) result(output)
  use string_module
  use err_module
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
  use err_module
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
end module

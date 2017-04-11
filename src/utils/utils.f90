! ======================================================================
! various utilities
! ======================================================================
module utils_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  interface format_path
    module procedure format_path_character
    module procedure format_path_String
  end interface
contains

! ----------------------------------------------------------------------
! Returns an array containing the command line arguments.
! ----------------------------------------------------------------------
function command_line_args() result(args)
  implicit none

  type(String), allocatable :: args(:)
  
  ! Temporary variables.
  integer         :: i,ialloc
  integer         :: arg_count
  character(1000) :: temp_char

  arg_count = command_argument_count()
  allocate (args(arg_count), stat=ialloc); call err(ialloc)
  do i=1,arg_count
    call get_command_argument(i, temp_char)
    args(i) = trim(temp_char)
  enddo
end function

! ----------------------------------------------------------------------
! Make a directory, if it doesn't already exist.
! ----------------------------------------------------------------------
subroutine mkdir(dirname)
  implicit none
  
  type(String), intent(in) :: dirname
  
  call system_call('if [ ! -e '//dirname//' ]; then mkdir '//dirname//'; fi')
end subroutine

! ----------------------------------------------------------------------
! Takes a directory name, and converts it into an absolute path
!    in standard format (without a trailing '/').
! ----------------------------------------------------------------------
function format_path_character(path,cwd) result(output)
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
! Vector L2 norm.
! Replicates norm2() from f2008 standard.
! ----------------------------------------------------------------------
pure function l2_norm(input) result(output)
  implicit none
  
  real(dp), intent(in) :: input(:)
  real(dp)             :: output
  
  output = dsqrt(dot_product(input,input))
end function

! ----------------------------------------------------------------------
! Vector outer product.
! ----------------------------------------------------------------------
function outer_product(input1,input2) result(output)
  implicit none
  
  real(dp), intent(in)  :: input1(:)
  real(dp), intent(in)  :: input2(:)
  real(dp), allocatable :: output(:,:)
  
  ! Temporary variables.
  integer :: i,ialloc
  
  allocate(output(size(input1), size(input2)), stat=ialloc); call err(ialloc)
  do i=1,size(input2)
    output(:,i) = input2(i) * input1
  enddo
end function
end module

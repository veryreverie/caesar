! ======================================================================
! various utilities
! ======================================================================
module utils_module
  use constants_module, only : dp
  use string_module
  use io_module
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
  allocate (args(arg_count+1), stat=ialloc); call err(ialloc)
  do i=0,arg_count
    call get_command_argument(i, temp_char)
    args(i+1) = trim(temp_char)
  enddo
end function

! ----------------------------------------------------------------------
! Make a directory, if it doesn't already exist.
! ----------------------------------------------------------------------
subroutine mkdir(dirname)
  implicit none
  
  type(String), intent(in) :: dirname
  
  integer :: result_code
  
  result_code = system_call( &
     & 'if [ ! -e '//dirname//' ]; then mkdir '//dirname//'; fi')
  if (result_code/=0) then
    call print_line('Error: failed to make directory: '//dirname)
    call err()
  endif
end subroutine

! ----------------------------------------------------------------------
! Converts a file seedname into the appropriate dft input or output filename.
! ----------------------------------------------------------------------
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
  
  output = sqrt(dot_product(input,input))
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

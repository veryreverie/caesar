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

! ----------------------------------------------------------------------
! Factorial.
! Stores calculated results to avoid excess computation.
! ----------------------------------------------------------------------
function factorial(input) result(output)
  implicit none
  
  integer, intent(in) :: input
  integer             :: output
  
  ! Lookup is saved between calls of this function.
  integer, allocatable, save :: lookup(:)
  integer, allocatable       :: temp(:)
  
  integer :: i,ialloc
  
  ! Initialise lookup the first time factorial is called.
  if (.not. allocated(lookup)) then
    lookup = [1]
  endif
  
  ! If the requested factorial has not previously been calculated,
  !    extends lookup to include it.
  if (input>=size(lookup)) then
    temp = lookup
    deallocate(lookup, stat=ialloc); call err(ialloc)
    allocate(lookup(input+1), stat=ialloc); call err(ialloc)
    lookup(:size(temp)) = temp
    do i=size(temp)+1,size(lookup)
      lookup(i) = (i-1)*lookup(i-1)
    enddo
  endif
  
  ! lookup(1) = 0!, so lookup(n+1) = n!.
  output = lookup(input+1)
end function
end module

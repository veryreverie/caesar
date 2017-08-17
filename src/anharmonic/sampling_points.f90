! ======================================================================
! Holds information about points to be sampled.
! ======================================================================
module sampling_points_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use coupling_module
  use linear_algebra_module
  
  ! The ids of the sampling point.
  ! The corresponding displacement depends on the chosen grid type.
  type :: SamplingPoint
    integer, allocatable :: indices(:)
    ! Whether or not this is a duplicate of another SamplingPoint, e.g.
    !    the point [17,0] in the coupling [1,2] is a duplicate of
    !    the point [17,0] in the coupling [1].
    logical              :: duplicate
  end type
  
  type :: CouplingSampling
    type(CoupledModes)               :: coupling
    type(SamplingPoint), allocatable :: sampling_points(:)
    real(dp),            allocatable :: energy(:)
    type(RealVector),    allocatable :: forces(:,:)
  end type
contains

subroutine write_sampling_points_file(this, filename)
  implicit none
  
  type(SamplingPoint), intent(in) :: this(:)
  type(String),        intent(in) :: filename
  
  integer :: sampling_file
  integer :: i
  
  sampling_file = open_write_file(filename)
  call print_line( sampling_file, &
                 & '! Sampling points in normal-mode co-ordinates.')
  call print_line( sampling_file, &
                 & '! Each co-ordinate is in multiples of delta.')
  do i=1,size(this)
    call print_line(sampling_file,this%duplicate//' '//this(i)%indices)
  enddo
  close(sampling_file)
end subroutine

function read_sampling_points_file(filename) result(this)
  implicit none
  
  type(String), intent(in)         :: filename
  type(SamplingPoint), allocatable :: this(:)
  
  type(String), allocatable :: sampling_file(:)
  
  integer                   :: i,ialloc
  type(String), allocatable :: line(:)
  
  sampling_file = read_lines(filename)
  
  allocate(this(size(sampling_file)-2), stat=ialloc); call err(ialloc)
  
  do i=3,size(sampling_file)
    line = split(sampling_file(i))
    this(i-2)%duplicate = lgcl(line(1))
    this(i-2)%indices = int(line(2:))
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the sampling points for a given coupling in the cubic grid.
! ----------------------------------------------------------------------
! The sampling points for coupling [1,3] has co-ordinates
!    [u1,0,u3,0,...,0]
!    for all u1 and u3 in set [-no_sampling_points,no_sampling_points]
! This function calls itself recursively, filling in indices from the right.
! In the [1,3] example, first the u3 indices are filled by a nested call,
!    and then the u1 indices are filled.
! Any sampling point with ui=0 for i in the coupling will be labled duplicate.
! If no_sampling_points=1, the result for coupling [1,3] will be
!
! [ [-1, 0, -1, 0, ..., 0] duplicate = .false. ,
!   [-1, 0,  0, 0, ..., 0] duplicate = .true.  ,
!   [-1, 0,  1, 0, ..., 0] duplicate = .false. ,
!   [ 0, 0, -1, 0, ..., 0] duplicate = .true.  ,
!   [ 0, 0,  0, 0, ..., 0] duplicate = .true.  ,
!   [ 0, 0,  1, 0, ..., 0] duplicate = .true.  ,
!   [ 1, 0, -1, 0, ..., 0] duplicate = .false. ,
!   [ 1, 0,  0, 0, ..., 0] duplicate = .true.  ,
!   [ 1, 0,  1, 0, ..., 0] duplicate = .false. ]
recursive function cubic_sampling_points(coupling, no_modes, &
   & no_sampling_points) result(output)
  use coupling_module
  use linear_algebra_module
  implicit none
  
  integer, intent(in)              :: coupling(:)
  integer, intent(in)              :: no_modes
  integer, intent(in)              :: no_sampling_points
  type(SamplingPoint), allocatable :: output(:)
  
  type(SamplingPoint), allocatable :: temp(:)
  
  integer :: i,i2,j,ialloc
  
  if (size(coupling)==0) then
    ! The base case: no modes, so only the point [0,0,0,...,0] is needed.
    allocate(output(1), stat=ialloc); call err(ialloc)
    output(1)%indices = int(zeroes(no_modes))
    output(1)%duplicate = .false.
  else
    ! Copy the sampling points corresponding to coupling(2:), [0,...,0,0,~~~]
    ! Make one copy for each i in set [-no_sampling_points,no_sampling_points]
    ! Set duplicate=.true. for i=0.
    ! If no_sampling_points=2 then the result is:
    !    [0,...,0,-2,~~~]
    !           ...
    !    [0,...,0,-2,~~~]
    !    [0,...,0,-1,~~~]
    !           ...
    !    [0,...,0,-1,~~~]
    !    [0,...,0, 0,~~~] (duplicate = .true.)
    !           ...
    !    [0,...,0, 0,~~~] (duplicate = .true.)
    !    [0,...,0, 1,~~~]
    !           ...
    !    [0,...,0, 1,~~~]
    !    [0,...,0, 2,~~~]
    !           ...
    !    [0,...,0, 2,~~~]
    temp = cubic_sampling_points(coupling(2:),no_modes,no_sampling_points)
    allocate(output(size(temp)*(2*no_sampling_points+1)), stat=ialloc)
    ! Loop from i2= -no_sampling_points to i2= no_sampling_points.
    do i=1,2*no_sampling_points+1
      i2 = i-no_sampling_points-1
      output((i-1)*size(temp)+1 : i*size(temp)) = temp
      do j=(i-1)*size(temp)+1,i*size(temp)
        output(j)%indices(coupling(1)) = i2
        if (i2==0) then
          output(j)%duplicate = .true.
        endif
      enddo
    enddo
  endif
end function
end module

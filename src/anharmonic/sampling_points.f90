! ======================================================================
! Holds information about points to be sampled.
! ======================================================================
module sampling_points_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! The ids of the sampling point.
  ! The corresponding displacement depends on the chosen grid type.
  type :: SamplingPoint
    integer, allocatable :: indices(:)
    ! Whether or not this is a duplicate of another SamplingPoint, e.g.
    !    the point [17,0] in the coupling [1,2] is a duplicate of
    !    the point [17,0] in the coupling [1].
    logical              :: duplicate
  end type
  
  ! Displacement vectors in normal mode co-ordinates.
  type :: DisplacementData
    real(dp), allocatable :: displacements(:)
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
  implicit none
  
  integer, intent(in)              :: coupling(:)
  integer, intent(in)              :: no_modes
  integer, intent(in)              :: no_sampling_points
  type(SamplingPoint), allocatable :: output(:)
  
  integer            :: no_points
  
  integer :: mid
  integer :: pos,neg
  integer :: i,j,ialloc
  
  allocate( output((2*no_sampling_points+1)**size(coupling)), &
          & stat=ialloc); call err(ialloc)
  
  if (size(coupling)==0) then
    ! The base case: no modes, so only the point [0,0,0,...,0] is needed.
    allocate(output(1)%indices(no_modes), stat=ialloc); call err(ialloc)
    output(1)%indices = 0
    output(1)%duplicate = .false.
  else
    ! Fill the middle with terms which look like (if coupling(1)=4)
    !    [0,0,0,0,[lower terms]]
    ! i.e. zeros up to and including the mode at coupling(1).
    no_points = (2*no_sampling_points+1)**(size(coupling(2:)))
    mid = no_sampling_points*no_points
    output(mid+1:mid+no_points) = cubic_sampling_points( coupling(2:), &
                                                       & no_modes,     &
                                                       & no_sampling_points)
    
    ! Copy the middle to each point, and then update the indices of
    !    mode coupling(1), to end up with e.g. if coupling(1)=4,
    !    [0,0,0,-2,[lower terms]]
    !    [0,0,0,-1,[lower terms]]
    !    [0,0,0, 0,[lower terms]] (duplicate = .true., set below)
    !    [0,0,0, 1,[lower terms]]
    !    [0,0,0, 2,[lower terms]]
    do i=1,no_sampling_points
      pos = mid + i*no_points
      neg = mid - i*no_points
      output(pos+1:pos+no_points) = output(mid+1:mid+no_points)
      output(neg+1:neg+no_points) = output(mid+1:mid+no_points)
      do j=1,no_points
        output(pos+j)%indices(coupling(1)) =  i
        output(neg+j)%indices(coupling(1)) = -i
      enddo
    enddo
    
    ! Set the duplicate=.true. on the points with indices(coupling(1))=0
    do i=1,no_points
      output(mid+i)%duplicate = .true.
    enddo
  endif
end function
end module

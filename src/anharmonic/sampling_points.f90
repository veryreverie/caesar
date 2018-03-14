! ======================================================================
! Holds information about points to be sampled.
! ======================================================================
module sampling_points_module
  use precision_module
  use io_module
  use algebra_module
  
  use coupling_module
  use normal_mode_module
  use mode_vector_module
  implicit none
  
  ! The ids of the sampling point.
  ! The corresponding displacement depends on the chosen grid type.
  type :: SamplingPoint
    integer, allocatable :: indices(:)
    type(RealModeVector) :: displacement
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
  
  type(OFile) :: sampling_file
  
  integer :: i
  
  sampling_file = filename
  call sampling_file%print_line( &
     & '! Sampling points in normal-mode co-ordinates.')
  do i=1,size(this)
    call sampling_file%print_line( '')
    call sampling_file%print_line( 'Indices      : '//this(i)%indices)
    call sampling_file%print_line( 'Displacement : '// &
                                 & this(i)%displacement%vector)
    call sampling_file%print_line( 'Duplicate    : '//this(i)%duplicate)
  enddo
end subroutine

function read_sampling_points_file(filename) result(this)
  implicit none
  
  type(String), intent(in)         :: filename
  type(SamplingPoint), allocatable :: this(:)
  
  type(IFile) :: sampling_file
  
  integer :: no_points
  
  integer                   :: i,ialloc
  type(String), allocatable :: line(:)
  
  sampling_file = filename
  
  no_points = (size(sampling_file)-1)/4
  allocate(this(no_points), stat=ialloc); call err(ialloc)
  
  do i=1,no_points
    line = split(sampling_file%line(4*(i-1)+3))
    this(i)%indices = int(line(3:))
    line = split(sampling_file%line(4*(i-1)+4))
    this(i)%displacement%vector = dble(line(3:))
    line = split(sampling_file%line(4*(i-1)+5))
    this(i)%duplicate = lgcl(line(3))
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the sampling points for a given coupling.
! ----------------------------------------------------------------------
function generate_sampling_points(grid_type,coupling,no_modes, &
   & no_sampling_points,sample_spacing) result(output)
  use coupling_module
  use grid_types_module
  implicit none
  
  type(String), intent(in)         :: grid_type
  integer,      intent(in)         :: coupling(:)
  integer,      intent(in)         :: no_modes
  integer,      intent(in)         :: no_sampling_points
  real(dp),     intent(in)         :: sample_spacing(:)
  type(SamplingPoint), allocatable :: output(:)
  
  integer :: i,j,ialloc
  
  integer, allocatable :: grid(:,:)
  
  ! Generate a grid of points corresponding to the sampling points.
  if (grid_type=='cubic') then
    grid = generate_cubic_grid(  size(coupling),     &
                              & -no_sampling_points, &
                              &  no_sampling_points)
  elseif (grid_type=='octahedral' .or. grid_type=='spherical') then
    grid = generate_octahedral_grid( size(coupling),     &
                                   & no_sampling_points, &
                                   & include_negatives=.true.)
  endif
  
  ! Generate output.
  allocate(output(size(grid,2)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%indices = int(zeroes(no_modes))
    output(i)%duplicate = .false.
    do j=1,size(coupling)
      output(i)%indices(coupling(j)) = grid(j,i)
      if (grid(j,i)==0) then
        output(i)%duplicate = .true.
      endif
    enddo 
    output(i)%displacement = calculate_displacement( grid_type,         &
                                                   & output(i)%indices, &
                                                   & sample_spacing)
  enddo
end function
end module

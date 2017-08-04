! ======================================================================
! Holds information about points to be sampled.
! ======================================================================
module sampling_points_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  type :: SamplingPoint
    integer, allocatable :: coordinates(:)
  end type
  
contains

subroutine write_sampling_points_file(this, filename)
  implicit none
  
  type(SamplingPoint), intent(in) :: this(:)
  type(String),        intent(in) :: filename
  
  integer :: sample_file
  integer :: i
  
  sample_file = open_write_file(filename)
  call print_line(sample_file,'! Sampling points in normal-mode co-ordinates.')
  call print_line(sample_file,'! Each co-ordinate is in multiples of delta.')
  do i=1,size(this)
    call print_line(sample_file,this(i)%coordinates)
  enddo
  close(sample_file)
end subroutine

function read_sampling_points_file(filename) result(this)
  implicit none
  
  type(String), intent(in)         :: filename
  type(SamplingPoint), allocatable :: this(:)
  
  type(String), allocatable :: sample_file(:)
  integer                   :: i,ialloc
  
  sample_file = read_lines(filename)
  
  allocate(this(size(sample_file)-2), stat=ialloc); call err(ialloc)
  
  do i=3,size(sample_file)
    this(i-2)%coordinates = int(split(sample_file(i)))
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the sampling points for a given coupling in the cubic grid.
! ----------------------------------------------------------------------
! The sampling points for [1,2] has co-ordinates(u1,u2,0,...,0), for all
!    u1 and u2 in set [-no_sampling_points,no_sampling_points] but /=0
recursive function cubic_sampling_points(coupling, no_modes, &
   & no_sampling_points) result(output)
  use coupling_module
  implicit none
  
  type(CoupledModes), intent(in)   :: coupling
  integer,            intent(in)   :: no_modes
  integer,            intent(in)   :: no_sampling_points
  type(SamplingPoint), allocatable :: output(:)
  
  integer            :: no_coupled
  type(CoupledModes) :: partial_coupling
  integer            :: no_points
  
  integer :: i,j,ialloc
  
  no_coupled = size(coupling%modes)
  allocate( output((2*no_sampling_points)**no_coupled), &
          & stat=ialloc); call err(ialloc)
  if (no_coupled==1) then
    do i=1,no_sampling_points
      allocate(output(i)%coordinates(no_modes), stat=ialloc); call err(ialloc)
      output(i)%coordinates = 0
      output(i)%coordinates(coupling%modes(1)) = -no_sampling_points-1+i
      
      allocate( output(no_sampling_points+i)%coordinates(no_modes), &
              & stat=ialloc); call err(ialloc)
      output(no_sampling_points+i)%coordinates = 0
      output(no_sampling_points+i)%coordinates(coupling%modes(1)) = i
    enddo
  else
    partial_coupling%modes = coupling%modes(2:)
    no_points = (2*no_sampling_points)**(no_coupled-1)
    output(:no_points) = cubic_sampling_points( partial_coupling, &
                                              & no_modes,         &
                                              & no_sampling_points)
    do i=2,no_sampling_points
      output(no_points*(i-1)+1:no_points*i) = output(:no_points)
    enddo
    do i=1,no_sampling_points
      do j=1,2*no_sampling_points
        output(no_points*(i-1)+j)%coordinates(coupling%modes(1)) = &
           & -no_sampling_points-1+i
        output(no_points*(i-1+no_sampling_points)+j)%coordinates( &
           & coupling%modes(1)) = i
      enddo
    enddo
  endif
end function
end module

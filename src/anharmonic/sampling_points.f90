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
  end type
  
  ! Displacement vectors in normal mode co-ordinates.
  type :: DisplacementData
    real(dp), allocatable :: displacements(:)
  end type
  
  interface cubic_sampling_points
    module procedure cubic_sampling_points_CoupledModes
    module procedure cubic_sampling_points_all_CoupledModes
  end interface
  
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
    call print_line(sampling_file,this(i)%indices)
  enddo
  close(sampling_file)
end subroutine

function read_sampling_points_file(filename) result(this)
  implicit none
  
  type(String), intent(in)         :: filename
  type(SamplingPoint), allocatable :: this(:)
  
  type(String), allocatable :: sampling_file(:)
  integer                   :: i,ialloc
  
  sampling_file = read_lines(filename)
  
  allocate(this(size(sampling_file)-2), stat=ialloc); call err(ialloc)
  
  do i=3,size(sampling_file)
    this(i-2)%indices = int(split(sampling_file(i)))
  enddo
end function

! ----------------------------------------------------------------------
! Calculates the sampling points for a given coupling in the cubic grid.
! ----------------------------------------------------------------------
! The sampling points for [1,2] has co-ordinates(u1,u2,0,...,0), for all
!    u1 and u2 in set [-no_sampling_points,no_sampling_points] but /=0
recursive function cubic_sampling_points_CoupledModes(coupling, no_modes, &
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
      allocate(output(i)%indices(no_modes), stat=ialloc); call err(ialloc)
      output(i)%indices = 0
      output(i)%indices(coupling%modes(1)) = -no_sampling_points-1+i
      
      allocate( output(no_sampling_points+i)%indices(no_modes), &
              & stat=ialloc); call err(ialloc)
      output(no_sampling_points+i)%indices = 0
      output(no_sampling_points+i)%indices(coupling%modes(1)) = i
    enddo
  else
    partial_coupling%modes = coupling%modes(2:)
    no_points = (2*no_sampling_points)**(no_coupled-1)
    output(:no_points) = cubic_sampling_points_CoupledModes( &
       & partial_coupling, &
       & no_modes,        &
       & no_sampling_points)
    do i=2,no_sampling_points
      output(no_points*(i-1)+1:no_points*i) = output(:no_points)
    enddo
    do i=1,no_sampling_points
      do j=1,2*no_sampling_points
        output(no_points*(i-1)+j)%indices(coupling%modes(1)) = &
           & -no_sampling_points-1+i
        output(no_points*(i-1+no_sampling_points)+j)%indices( &
           & coupling%modes(1)) = i
      enddo
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Calculates the sampling points for all couplings in the cubic grid.
! ----------------------------------------------------------------------
recursive function cubic_sampling_points_all_CoupledModes(coupling, no_modes, &
   & no_sampling_points) result(output)
  use coupling_module
  implicit none
  
  type(CoupledModes), intent(in)   :: coupling(:)
  integer,            intent(in)   :: no_modes
  integer,            intent(in)   :: no_sampling_points
  type(SamplingPoint), allocatable :: output(:)
  
  integer, allocatable :: points_per_coupling(:)
  
  integer :: i,j,ialloc
  
  ! Calculate the number of sampling points associated with each coupling.
  allocate( points_per_coupling(size(coupling)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(coupling)
    points_per_coupling(i) = (2*no_sampling_points)**size(coupling(i)%modes)
  enddo
  
  ! Calculate sampling points.
  
  ! The first point is the unperturbed supercell.
  allocate( output(sum(points_per_coupling)+1), &
          & stat=ialloc); call err(ialloc)
  allocate(output(1)%indices(no_modes), stat=ialloc); call err(ialloc)
  output(1)%indices = 0
  
  ! All further points are organised in order of coupling.
  j = 1
  do i=1,size(coupling)
    output(j+1:j+points_per_coupling(i)) =   &
       & cubic_sampling_points( coupling(i), &
       &                        no_modes,    &
       &                        no_sampling_points)
    j = j + points_per_coupling(i)
  enddo
end function
end module

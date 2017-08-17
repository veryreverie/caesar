! ======================================================================
! Anharmonic calculations for different grid types.
! ======================================================================
module grid_types_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  private
  
  public :: calculate_displacement
  public :: calculate_cartesian_displacement
contains

! ----------------------------------------------------------------------
! Calculates the displacement associated with a sampling point.
! Returns the output in normal mode co-ordinates.
! ----------------------------------------------------------------------
function calculate_displacement(grid_type,sampling_point,sample_spacing) &
   & result(output)
  use sampling_points_module
  use normal_mode_module
  type(String),        intent(in) :: grid_type
  type(SamplingPoint), intent(in) :: sampling_point
  real(dp),            intent(in) :: sample_spacing(:)
  type(ModeVector)                :: output
  
  ! Calculate displacement in normal mode co-ordinates,
  !    depending on grid type.
  if (grid_type=='cubic') then
    output = calculate_displacement_cubic(sampling_point, sample_spacing)
  elseif (grid_type=='octahedral') then
    output = calculate_displacement_octahedral(sampling_point, sample_spacing)
  elseif (grid_type=='spherical') then
    output = calculate_displacement_spherical(sampling_point, sample_spacing)
  else
    call err()
  endif
end function

function calculate_displacement_cubic(sampling_point,sample_spacing) &
   & result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  type(SamplingPoint), intent(in) :: sampling_point
  real(dp),            intent(in) :: sample_spacing(:)
  type(ModeVector)                :: output
  
  ! Temporary variables.
  integer :: i,ialloc
  
  allocate( output%vector(size(sampling_point%indices)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(output%vector)
    output%vector(i) = sampling_point%indices(i) * sample_spacing(i)
  enddo
end function

function calculate_displacement_octahedral(sampling_point,sample_spacing) &
   & result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  type(SamplingPoint), intent(in) :: sampling_point
  real(dp),            intent(in) :: sample_spacing(:)
  type(ModeVector)                :: output
  
  ! The octahedral and cubic sampling points are indexed in the same way.
  output = calculate_displacement_cubic(sampling_point,sample_spacing)
end function

function calculate_displacement_spherical(sampling_point,sample_spacing) &
   & result(output)
  use sampling_points_module
  use normal_mode_module
  implicit none
  
  type(SamplingPoint), intent(in) :: sampling_point
  real(dp),            intent(in) :: sample_spacing(:)
  type(ModeVector)                :: output
  
  real(dp) :: radius
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Calculate radius.
  radius = sum(sampling_point%indices)
  
  allocate( output%vector(size(sampling_point%indices)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(output%vector)
    output%vector(i) = sin(sampling_point%indices(i)/radius) &
                   & * radius                                &
                   & * sample_spacing(i)
  enddo
end function
end module

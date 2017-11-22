! ======================================================================
! Anharmonic calculations for different grid types.
! ======================================================================
module grid_types_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  private
  
  public :: calculate_displacement
  public :: generate_cubic_grid
  public :: generate_octahedral_grid
  public :: octahedral_grid_size
  
  type, public :: PointData
    integer, allocatable :: point(:)
    logical              :: duplicate
  end type
contains

! ----------------------------------------------------------------------
! Calculates the displacement associated with a sampling point.
! Returns the output in normal mode co-ordinates.
! ----------------------------------------------------------------------
function calculate_displacement(grid_type,sampling_point_indices, &
   & sample_spacing,qpoint) result(output)
  use normal_mode_module
  use qpoints_module
  type(String),     intent(in) :: grid_type
  integer,          intent(in) :: sampling_point_indices(:)
  real(dp),         intent(in) :: sample_spacing(:)
  type(QpointData), intent(in) :: qpoint
  type(RealModeVector)         :: output
  
  ! Calculate displacement in normal mode co-ordinates,
  !    depending on grid type.
  if (grid_type=='cubic') then
    output = calculate_displacement_cubic( sampling_point_indices, &
                                         & sample_spacing,         &
                                         & qpoint)
  elseif (grid_type=='octahedral') then
    output = calculate_displacement_octahedral( sampling_point_indices, &
                                              & sample_spacing,         &
                                              & qpoint)
  elseif (grid_type=='spherical') then
    output = calculate_displacement_spherical( sampling_point_indices, &
                                             & sample_spacing,         &
                                             & qpoint)
  else
    call err()
  endif
end function

function calculate_displacement_cubic(sampling_point_indices, &
   & sample_spacing,qpoint) result(output)
  use normal_mode_module
  use qpoints_module
  implicit none
  
  integer,          intent(in) :: sampling_point_indices(:)
  real(dp),         intent(in) :: sample_spacing(:)
  type(QpointData), intent(in) :: qpoint
  type(RealModeVector) :: output
  
  integer              :: no_modes
  
  ! Temporary variables.
  integer :: i,j,mode,ialloc
  
  no_modes = size(sampling_point_indices)
  
  if (qpoint%is_paired_qpoint) then
    allocate(output%vector(no_modes,1), stat=ialloc); call err(ialloc)
    do mode=1,no_modes
      output%vector(mode,1) = sampling_point_indices(mode) &
                          & * sample_spacing(mode)
    enddo
  else
    allocate(output%vector(no_modes/2,2), stat=ialloc); call err(ialloc)
    do i=1,no_modes/2
      do j=1,2
        mode = j*no_modes/2+i-1
        output%vector(i,j) = sampling_point_indices(mode) &
                         & * sample_spacing(mode)
      enddo
    enddo
  endif
end function

function calculate_displacement_octahedral(sampling_point_indices, &
   & sample_spacing,qpoint) result(output)
  use normal_mode_module
  use qpoints_module
  implicit none
  
  integer,          intent(in) :: sampling_point_indices(:)
  real(dp),         intent(in) :: sample_spacing(:)
  type(QpointData), intent(in) :: qpoint
  type(RealModeVector)         :: output
  
  ! The octahedral and cubic sampling points are indexed in the same way.
  output = calculate_displacement_cubic( sampling_point_indices, &
                                       & sample_spacing,         &
                                       & qpoint)
end function

function calculate_displacement_spherical(sampling_point_indices, &
   & sample_spacing,qpoint) result(output)
  use normal_mode_module
  use qpoints_module
  implicit none
  
  integer,          intent(in) :: sampling_point_indices(:)
  real(dp),         intent(in) :: sample_spacing(:)
  type(QpointData), intent(in) :: qpoint
  type(RealModeVector)         :: output
  
  integer :: radius
  integer :: no_modes
  
  ! Temporary variables.
  integer :: i,j,mode,ialloc
  
  no_modes = size(sampling_point_indices)
  
  ! Calculate radius.
  radius = sum(sampling_point_indices)
  
  if (qpoint%is_paired_qpoint) then
    allocate(output%vector(no_modes,1), stat=ialloc); call err(ialloc)
    if (radius==0) then
      output%vector = 0
    else
      do mode=1,no_modes
        output%vector(mode,1) = sin( sampling_point_indices(mode) &
                          &        / real(radius,dp) )            &
                          & * radius                              &
                          & * sample_spacing(mode)
      enddo
    endif
  else
    allocate(output%vector(no_modes/2,2), stat=ialloc); call err(ialloc)
    if (radius==0) then
      output%vector = 0
    else
      do i=1,no_modes/2
        do j=1,2
          mode = j*no_modes/2+i-1
          output%vector(i,j) = sin( sampling_point_indices(mode) &
                            &        / real(radius,dp) )            &
                            & * radius                              &
                            & * sample_spacing(mode)
        enddo
      enddo
    endif
  endif
end function

! ----------------------------------------------------------------------
! Generates a grid of points.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Generates all points in a (hyper-)cubic grid.
! --------------------------------------------------
! Each point has no_dimensions elements.
! Each element is in the set [lower_bound,upper_bound].
! e.g. if no_dimensions=3, lower_bound=-1 and upper_bound=2 then output=
!
! [ [-1,-1,-1],
!   [-1,-1, 0],
!   [-1,-1, 1],
!   [-1,-1, 2],
!   [-1, 0,-1],
!       ...
!   [-1, 0, 2],
!   [-1, 1,-1],
!       ...
!   [-1, 2, 2],
!   [ 0,-1,-1]
!       ...
!   [ 2, 2, 2] ]
recursive function generate_cubic_grid(no_dimensions,lower_bound,upper_bound) &
   & result(output)
  implicit none
  
  integer, intent(in)  :: no_dimensions
  integer, intent(in)  :: lower_bound
  integer, intent(in)  :: upper_bound
  integer, allocatable :: output(:,:)
  
  ! The number of points between lower_bound and upper_bound inclusive.
  integer :: range
  
  ! The result of the recursive call, and its size.
  integer, allocatable :: smaller_grid(:,:)
  integer              :: small
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  range = upper_bound-lower_bound+1 ! No. point
  
  allocate( output(no_dimensions, range**no_dimensions), &
          & stat=ialloc); call err(ialloc)
  
  if (no_dimensions==0) then
    output = 0 ! This will be a 0 by 1 array.
  else
    smaller_grid = generate_cubic_grid(no_dimensions-1,lower_bound,upper_bound)
    small = size(smaller_grid,2)
    j = 0
    do i=lower_bound,upper_bound
      output(1,  j+1:j+small) = i
      output(2:, j+1:j+small) = smaller_grid
      j = j+small
    enddo
  endif
end function

! --------------------------------------------------
! Generates all points in a (hyper-)octahedral grid.
! --------------------------------------------------
! As generate_cubic_grid, except only points which satisfy
!    sum(abs(output(:,i))) <= upper_bound are included.
!
! lower_bound from generate_cubic_grid is replaced by include_negatives.
!    - if include_negatives=.false. then
!       lower_bound=0
!
!                 (no_dimensions  + upper_bound )!
!       There are -------------------------------- points.
!                 (no_dimensions! * upper_bound!)
!
!    - if include_negatives=.true. then
!       lower_bound=-upper_bound
!
!                 min(u_b,n_d)        u_b!         *  n_d!
!       There are    sum       2^i * -----------------------------
!                    i=0             (u_b-i)! * i! * (n_d-i)! * i!
recursive function generate_octahedral_grid(no_dimensions,upper_bound, &
   & include_negatives) result(output)
  implicit none
  
  integer, intent(in)  :: no_dimensions
  integer, intent(in)  :: upper_bound
  logical, intent(in)  :: include_negatives
  integer, allocatable :: output(:,:)
  
  ! The result of the recursive call, and its size.
  integer, allocatable :: smaller_grid(:,:)
  integer              :: small
  
  integer :: no_points
  integer :: min_i
  
  integer :: i,j,ialloc
  
  if (no_dimensions==0 .or. upper_bound==0) then
    allocate(output(no_dimensions,1), stat=ialloc); call err(ialloc)
    output = 0
  else
    no_points = octahedral_grid_size( no_dimensions, &
                                    & upper_bound,   &
                                    & include_negatives)
    allocate(output(no_dimensions,no_points), stat=ialloc); call err(ialloc)
    
    if (include_negatives) then
      min_i = -upper_bound
    else
      min_i = 0
    endif
    
    j = 0
    do i=min_i,upper_bound
      smaller_grid = generate_octahedral_grid( no_dimensions-1,    &
                                             & upper_bound-abs(i), &
                                             & include_negatives)
      small = size(smaller_grid,2)
      
      output(1,  j+1:j+small) = i
      output(2:, j+1:j+small) = smaller_grid
      j = j+small
    enddo
    
    ! Check that the output has been filled correctly.
    if (j/=no_points) then
      call err()
    endif
  endif
end function

! --------------------------------------------------
! Returns the number of points in an octahedral grid with given parameters.
! --------------------------------------------------
function octahedral_grid_size(no_dimensions,upper_bound,include_negatives) &
   & result(output)
  use utils_module, only : factorial
  implicit none
  
  integer, intent(in) :: no_dimensions
  integer, intent(in) :: upper_bound
  logical, intent(in) :: include_negatives
  integer             :: output
  
  integer :: i
  
  if (include_negatives) then
    output = 0
    do i=0,min(no_dimensions,upper_bound)
      output = output                       &
           & + 2**i                         &
           & * factorial(no_dimensions)     &
           & * factorial(upper_bound)       &
           & / ( factorial(i)**2            &
           &   * factorial(no_dimensions-i) &
           &   * factorial(upper_bound-i) )
    enddo
  else
    output = factorial(no_dimensions+upper_bound) &
         & / ( factorial(no_dimensions)           &
         &   * factorial(upper_bound) )
  endif
end function
end module

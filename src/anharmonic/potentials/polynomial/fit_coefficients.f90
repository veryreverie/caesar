! ======================================================================
! Routines for fitting the coefficients of polynomial basis functions.
! ======================================================================
module fit_coefficients_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
  use sampling_points_module
  use sample_result_module
  implicit none
  
  private
  
  public :: fit_coefficients
contains

! Uses L2 regression to calculate the coefficients of a set of basis functions.
function fit_coefficients(basis_functions,sampling_points,sample_results, &
   & modes,energy_force_ratio,potential) result(output)
  implicit none
  
  type(BasisFunction),        intent(in)           :: basis_functions(:)
  type(RealModeDisplacement), intent(in)           :: sampling_points(:)
  type(SampleResult),         intent(in)           :: sample_results(:)
  type(RealMode),             intent(in)           :: modes(:)
  real(dp),                   intent(in)           :: energy_force_ratio
  class(PotentialData),       intent(in), optional :: potential
  real(dp), allocatable                            :: output(:)
  
  real(dp)            :: energy
  type(RealModeForce) :: forces
  
  integer               :: dimensions
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  integer :: i,j,ialloc
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Each calculation yields size(modes) forces and one energy.
  dimensions = 1+size(modes)
  
  ! Calculate the energies and forces due to each basis function at each
  !    sampling point.
  a = construct_sample_matrix( basis_functions,   &
                             & sampling_points,   &
                             & modes,             &
                             & energy_force_ratio )
  
  ! Calculate the energies and forces sampled at each sampling point.
  b = construct_sample_vector( sampling_points,   &
                             & sample_results,    &
                             & potential,         &
                             & modes,             &
                             & energy_force_ratio )
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  output = dble(linear_least_squares(a, b))
end function

! Convert an energy and force into a single vector which can be inserted
!    into a matrix for passing to LAPACK.
function make_vector(energy,force,modes,energy_force_ratio) result(output)
  implicit none
  
  real(dp),            intent(in) :: energy
  type(RealModeForce), intent(in) :: force
  type(RealMode),      intent(in) :: modes(:)
  real(dp),            intent(in) :: energy_force_ratio
  real(dp), allocatable           :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(1+size(modes)), stat=ialloc); call err(ialloc)
  output(1) = energy / energy_force_ratio
  do i=1,size(modes)
    output(1+i) = force%force(modes(i))
  enddo
end function
end module

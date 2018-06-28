! ======================================================================
! Routines for fitting the coefficients of polynomial basis functions.
! ======================================================================
module fit_coefficients_module
  use common_module
  
  use anharmonic_common_module
  
  use basis_function_module
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
  
  real(dp),            allocatable :: basis_function_energies(:,:)
  type(RealModeForce), allocatable :: basis_function_forces(:,:)
  real(dp),            allocatable :: sample_energies(:)
  type(RealModeForce), allocatable :: sample_forces(:)
  
  integer               :: degrees_of_freedom
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  integer :: i,j,k,ialloc
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Calculate the energies and forces due to each basis function at each
  !    sampling point.
  allocate( basis_function_energies( size(sampling_points),  &
          &                          size(basis_functions)), &
          & basis_function_forces( size(sampling_points),    &
          &                        size(basis_functions)),   &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    do j=1,size(sampling_points)
      basis_function_energies(j,i) = basis_functions(i)%energy( &
                                           & sampling_points(j) )
      basis_function_forces(j,i) = basis_functions(i)%force(sampling_points(j))
    enddo
  enddo
  
  ! Calculate the energies and forces sampled at each sampling point.
  allocate( sample_energies(size(sampling_points)), &
          & sample_forces(size(sampling_points)),   &
          & stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    sample_energies(i) = sample_results(i)%energy
    sample_forces(i) = sample_results(i)%force
  enddo
  
  ! Subtract the energy and forces from the existing potential.
  if (present(potential)) then
    do i=1,size(sampling_points)
      sample_energies(i) = sample_energies(i) &
                       & - potential%energy(sampling_points(i))
      sample_forces(i) = sample_forces(i) &
                     & - potential%force(sampling_points(i))
    enddo
  endif
  
  ! Convert basis function energies and forces into a single matrix,
  !    and sampled energies and forces into a single vector.
  degrees_of_freedom = size(sampling_points) * (1+size(modes))
  allocate( a(degrees_of_freedom,size(basis_functions)), &
          & b(degrees_of_freedom),                       &
          & stat=ialloc); call err(ialloc)
  
  do i=1,size(basis_functions)
    k = 1
    do j=1,size(sampling_points)
      a(k:k+size(modes),i) = make_vector( basis_function_energies(j,i), &
                                        & basis_function_forces(j,i),   &
                                        & modes,                        &
                                        & energy_force_ratio)
      k = k + 1+size(modes)
    enddo
  enddo
  
  k = 1
  do i=1,size(sampling_points)
    b(k:k+size(modes)) = make_vector( sample_energies(i), &
                                    & sample_forces(i),   &
                                    & modes,              &
                                    & energy_force_ratio)
    k = k + 1+size(modes)
  enddo
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  output = dble(linear_least_squares(a,b))
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
  
  integer :: i,j,ialloc
  
  allocate(output(1+size(modes)), stat=ialloc); call err(ialloc)
  output(1) = energy / energy_force_ratio
  do i=1,size(modes)
    j = first(force%vectors%id==modes(i)%id, default=0)
    if (j==0) then
      output(1+i) = 0.0_dp
    else
      output(1+i) = force%vectors(j)%magnitude
    endif
  enddo
end function
end module

! ======================================================================
! Fits the coefficients of polynomial stress basis functions.
! ======================================================================
module fit_stress_coefficients_module
  use common_module
  
  use anharmonic_common_module
  
  use stress_basis_function_module
  use sample_result_module
  implicit none
  
  private
  
  public :: fit_stress_coefficients
contains

! Uses L2 regression to calculate the coefficients of a set of stress basis
!    functions.
function fit_stress_coefficients(basis_functions,sampling_points, &
   & sample_results,stress) result(output)
  implicit none
  
  type(StressBasisFunction),  intent(in)           :: basis_functions(:)
  type(RealModeDisplacement), intent(in)           :: sampling_points(:)
  type(SampleResult),         intent(in)           :: sample_results(:)
  class(StressData),          intent(in), optional :: stress
  real(dp), allocatable                            :: output(:)
  
  type(RealMatrix) :: stress_tensor
  
  real(dp), allocatable :: a(:,:)
  real(dp), allocatable :: b(:)
  
  integer :: i,j,ialloc
  
  ! Check inputs are consistent.
  if (size(sampling_points)/=size(sample_results)) then
    call print_line(CODE_ERROR//': The number of sampling points does not &
       &match the number of results.')
    call err()
  endif
  
  ! Calculate the stress due to each basis function at each sampling point.
  allocate( a(9*size(sampling_points), size(basis_functions)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(basis_functions)
    do j=1,size(sampling_points)
      stress_tensor = basis_functions(i)%stress(sampling_points(j))
      a((j-1)*9+1:j*9, i) = reshape(dble(stress_tensor), [9])
    enddo
  enddo
  
  ! Calculate the stress sampled at each sampling point.
  allocate(b(9*size(sampling_points)), stat=ialloc); call err(ialloc)
  do i=1,size(sampling_points)
    stress_tensor = sample_results(i)%stress()
    
    ! Subtract the existing stress.
    if (present(stress)) then
      stress_tensor = stress_tensor - stress%stress(sampling_points(i))
    endif
    
    b((i-1)*9+1:i*9) = reshape(dble(stress_tensor), [9])
  enddo
  
  ! Run linear least squares to get the basis function coefficients.
  ! This finds x s.t. (a.x-b)^2 is minimised.
  output = dble(linear_least_squares(a, b))
end function
end module

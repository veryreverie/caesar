! ======================================================================
! Various matrix tests.
! ======================================================================
module tests_submodule
  use precision_module
  use io_module
  
  use linear_algebra_submodule
  use algebra_utils_submodule
  implicit none
  
  private
  
  public :: check_symmetric
  public :: check_hermitian
  public :: check_orthogonal
  public :: check_unitary
  
  public :: check_orthonormal
  
  interface check_orthonormal
    module procedure check_orthonormal_RealMatrix
    module procedure check_orthonormal_ComplexMatrix
  end interface
contains

subroutine check_symmetric(input,logfile,warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)           :: input
  type(OFile),      intent(inout)        :: logfile
  real(dp),         intent(in), optional :: warning_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  real(dp) :: scale
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-10_dp
  endif
  
  scale = sum_squares(input)
  error = sum_squares(input-transpose(input))
  if (scale<1e-200_dp) then
    ! If scale is too small, error/scale will overflow.
    call logfile%print_line('||M|| too small to calculate fractional error in &
       &symmetry.')
  else
    error = sqrt(error/scale)
    call logfile%print_line('Fractional L2 error in symmetry: '//error)
    if (error>threshold) then
      call print_line(WARNING//': Matrix not as symmetric as expected. &
         &Fractional L2 error in symmetry: '//error)
    endif
  endif
end subroutine

subroutine check_hermitian(input,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)           :: input
  type(OFile),         intent(inout)        :: logfile
  real(dp),            intent(in), optional :: warning_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  real(dp) :: scale
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-10_dp
  endif
  
  scale = sum_squares(input)
  error = sum_squares(input-hermitian(input))
  if (scale<1e-200_dp) then
    ! If scale is too small, error/scale will overflow.
    call logfile%print_line('||M|| too small to calculate fractional error in &
       &Hermicity.')
  else
    error = sqrt(error/scale)
    call logfile%print_line('Fractional L2 error in Hermicity: '//error)
    if (error>threshold) then
      call print_line(WARNING//': Matrix not as Hermitian as expected. &
         &Fractional L2 error in Hermicity: '//error)
    endif
  endif
end subroutine

subroutine check_orthogonal(input,logfile,warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)           :: input
  type(OFile),      intent(inout)        :: logfile
  real(dp),         intent(in), optional :: warning_threshold
  
  real(dp) :: threshold
  
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-10_dp
  endif
  
  error = sum_squares( input*transpose(input) &
                   & - make_identity_matrix(size(input,1)))
  error = sqrt(error)
  call logfile%print_line('Fractional L2 error in orthogonality: '//error)
  if (error>threshold) then
    call print_line(WARNING//': Matrix not as orthogonal as expected. &
       &Fractional L2 error in orthogonality: '//error)
  endif
end subroutine

subroutine check_unitary(input,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)           :: input
  type(OFile),         intent(inout)        :: logfile
  real(dp),            intent(in), optional :: warning_threshold
  
  real(dp) :: threshold
  
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-10_dp
  endif
  
  error = sum_squares( input*hermitian(input) &
                   & - cmplxmat(make_identity_matrix(size(input,1))))
  error = sqrt(error)
  call logfile%print_line('Fractional L2 error in unitarity: '//error)
  if (error>threshold) then
    call print_line(WARNING//': Matrix not as unitary as expected. &
       &Fractional L2 error in unitarity: '//error)
  endif
end subroutine

subroutine check_orthonormal_RealMatrix(input,logfile,warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)           :: input
  type(OFile),      intent(inout)        :: logfile
  real(dp),         intent(in), optional :: warning_threshold
  
  real(dp) :: threshold
  
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-10_dp
  endif
  
  if (size(input,1)<size(input,2)) then
    error = sum_squares( input*transpose(input) &
                     & - cmplxmat(make_identity_matrix(size(input,1))))
  else
    error = sum_squares( transpose(input)*input &
                     & - cmplxmat(make_identity_matrix(size(input,2))))
  endif
  error = sqrt(error)
  call logfile%print_line('Fractional L2 error in unitarity: '//error)
  if (error>threshold) then
    call print_line(WARNING//': Matrix not as unitary as expected. &
       &Fractional L2 error in unitarity: '//error)
  endif
end subroutine

subroutine check_orthonormal_ComplexMatrix(input,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)           :: input
  type(OFile),         intent(inout)        :: logfile
  real(dp),            intent(in), optional :: warning_threshold
  
  real(dp) :: threshold
  
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-10_dp
  endif
  
  if (size(input,1)<size(input,2)) then
    error = sum_squares( input*hermitian(input) &
                     & - cmplxmat(make_identity_matrix(size(input,1))))
  else
    error = sum_squares( hermitian(input)*input &
                     & - cmplxmat(make_identity_matrix(size(input,2))))
  endif
  error = sqrt(error)
  call logfile%print_line('Fractional L2 error in unitarity: '//error)
  if (error>threshold) then
    call print_line(WARNING//': Matrix not as unitary as expected. &
       &Fractional L2 error in unitarity: '//error)
  endif
end subroutine
end module

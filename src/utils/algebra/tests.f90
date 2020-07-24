! ======================================================================
! Various matrix tests.
! ======================================================================
module tests_module
  use precision_module
  use io_module
  
  use linear_algebra_module
  use algebra_utils_module
  implicit none
  
  private
  
  public :: check_identity
  public :: check_symmetric
  public :: check_hermitian
  public :: check_orthogonal
  public :: check_unitary
  public :: check_orthonormal
  public :: check_real
  public :: check_int
  
  interface check_identity
    module procedure check_identity_String
    module procedure check_identity_character
  end interface
  
  interface check_symmetric
    module procedure check_symmetric_String
    module procedure check_symmetric_character
  end interface
  
  interface check_hermitian
    module procedure check_hermitian_String
    module procedure check_hermitian_character
  end interface
  
  interface check_orthogonal
    module procedure check_orthogonal_String
    module procedure check_orthogonal_character
  end interface
  
  interface check_unitary
    module procedure check_unitary_String
    module procedure check_unitary_character
  end interface
  
  interface check_orthonormal
    module procedure check_orthonormal_RealMatrix_String
    module procedure check_orthonormal_RealMatrix_character
    module procedure check_orthonormal_ComplexMatrix_String
    module procedure check_orthonormal_ComplexMatrix_character
  end interface
  
  interface check_real
    module procedure check_real_String
    module procedure check_real_character
  end interface
  
  interface check_int
    module procedure check_int_String
    module procedure check_int_character
  end interface
contains

! ----------------------------------------------------------------------
! Various checks for matrices.
! ----------------------------------------------------------------------
! Checks will print a warning message if they fail.
! matrix_name will be used to make messages more descriptive.
! If a logfile is present then it will always be written to.
! warning_threshold is the threshold at which the test is failed.
! ignore_threshold is the norm of the matrix below which the test is ignored.

subroutine check_identity_String(input,matrix_name,logfile, &
   & warning_threshold,ignore_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  type(String),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  real(dp),         intent(in),    optional :: ignore_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  
  if (size(input,1)/=size(input,2)) then
    call print_line(ERROR//': Matrix cannot be identity because it is not &
       &square.')
    call err()
  endif
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  error = sqrt(sum_squares(input-make_identity_matrix(size(input,1))))
  if (present(logfile)) then
    call logfile%print_line('L2 difference between identity and '// &
       & matrix_name//': '//error)
  endif
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' is not the identity &
       &to within the expected tolerance. L2 error : '//error)
  endif
end subroutine

subroutine check_identity_character(input,matrix_name,logfile, &
   & warning_threshold,ignore_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  character(*),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  real(dp),         intent(in),    optional :: ignore_threshold
  
  call check_identity( input,             &
                     & str(matrix_name),  &
                     & logfile,           &
                     & warning_threshold, &
                     & ignore_threshold)
end subroutine

subroutine check_symmetric_String(input,matrix_name,logfile, &
   & warning_threshold,ignore_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  type(String),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  real(dp),         intent(in),    optional :: ignore_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  real(dp) :: scale
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  scale = sqrt(sum_squares(input))
  
  ! Check if the norm of the input is smaller than ignore_threshold.
  if (present(ignore_threshold)) then
    if (scale<ignore_threshold) then
      if (present(logfile)) then
        call logfile%print_line(matrix_name//' has norm smaller than &
           &threshold for ignoring symmetry check.')
      endif
      return
    endif
  endif
  
  if (scale<1e-200_dp) then
    ! If scale is too small, error/scale will overflow.
    if (present(logfile)) then
      call logfile%print_line(matrix_name//' has norm too small to calculate &
         &fractional error in symmetry.')
    endif
  else
    error = sqrt(sum_squares(input-transpose(input))) / scale
    if (present(logfile)) then
      call logfile%print_line('Fractional L2 error in symmetry of '// &
         & matrix_name//': '//error)
    endif
    if (error>threshold) then
      call print_line(WARNING//': '//matrix_name//' not as symmetric as &
         &expected. Fractional L2 error in symmetry: '//error)
    endif
  endif
end subroutine

subroutine check_symmetric_character(input,matrix_name,logfile, &
   & warning_threshold,ignore_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  character(*),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  real(dp),         intent(in),    optional :: ignore_threshold
  
  call check_symmetric( input,             &
                      & str(matrix_name),  &
                      & logfile,           &
                      & warning_threshold, &
                      & ignore_threshold)
end subroutine

subroutine check_hermitian_String(input,matrix_name,logfile, &
   & warning_threshold,ignore_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  type(String),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  real(dp),            intent(in),    optional :: ignore_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  real(dp) :: scale
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  scale = sqrt(sum_squares(input))
  
  ! Check if the norm of the input is smaller than ignore_threshold.
  if (present(ignore_threshold)) then
    if (scale<ignore_threshold) then
      if (present(logfile)) then
        call logfile%print_line(matrix_name//' has norm smaller than &
           &threshold for ignoring symmetry check.')
      endif
      return
    endif
  endif
  
  if (scale<1e-200_dp) then
    ! If scale is too small, error/scale will overflow.
    if (present(logfile)) then
      call logfile%print_line(matrix_name//' has norm too small to calculate &
         &fractional error in symmetry.')
    endif
  else
    error = sqrt(sum_squares(input-hermitian(input))) / scale
    if (present(logfile)) then
      call logfile%print_line('Fractional L2 error in symmetry of '// &
         & matrix_name//': '//error)
    endif
    if (error>threshold) then
      call print_line(WARNING//': '//matrix_name//' not as symmetric as &
         &expected. Fractional L2 error in symmetry: '//error)
    endif
  endif
end subroutine

subroutine check_hermitian_character(input,matrix_name,logfile, &
   & warning_threshold,ignore_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  character(*),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  real(dp),            intent(in),    optional :: ignore_threshold
  
  call check_hermitian( input,             &
                      & str(matrix_name),  &
                      & logfile,           &
                      & warning_threshold, &
                      & ignore_threshold)
end subroutine

subroutine check_orthogonal_String(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  type(String),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  
  real(dp) :: threshold
  
  real(dp) :: error
  
  if (size(input,1)/=size(input,2)) then
    call print_line(ERROR//': '//matrix_name//" cannot be orthogonal because &
       &it isn't square.")
    call err()
  endif
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  error = sqrt(sum_squares( input*transpose(input) &
                        & - make_identity_matrix(size(input,1))))
  
  if (present(logfile)) then
    call logfile%print_line('Fractional L2 error in orthogonality of '// &
       & matrix_name//': '//error)
  endif
  
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' not as orthogonal as &
       &expected. Fractional L2 error in orthogonality: '//error)
  endif
end subroutine

subroutine check_orthogonal_character(input,matrix_name,logfile, &
   & warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  character(*),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  
  call check_orthogonal(input,str(matrix_name),logfile,warning_threshold)
end subroutine

subroutine check_unitary_String(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  type(String),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  
  real(dp) :: threshold
  
  real(dp) :: error
  
  if (size(input,1)/=size(input,2)) then
    call print_line(ERROR//': '//matrix_name//" cannot be unitary because it &
       &isn't square.")
    call err()
  endif
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  error = sqrt(sum_squares( input*hermitian(input) &
                        & - cmplxmat(make_identity_matrix(size(input,1)))))
  
  if (present(logfile)) then
    call logfile%print_line('Fractional L2 error in unitarity of '// &
       & matrix_name//': '//error)
  endif
  
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' not as unitary as expected. &
       &Fractional L2 error in unitarity: '//error)
  endif
end subroutine

subroutine check_unitary_character(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  character(*),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  
  call check_unitary(input,str(matrix_name),logfile,warning_threshold)
end subroutine

subroutine check_orthonormal_RealMatrix_String(input,matrix_name,logfile, &
   & warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  type(String),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  if (size(input,1)<size(input,2)) then
    error = sum_squares( input*transpose(input) &
                     & - cmplxmat(make_identity_matrix(size(input,1))))
  else
    error = sum_squares( transpose(input)*input &
                     & - cmplxmat(make_identity_matrix(size(input,2))))
  endif
  error = sqrt(error)
  
  if (present(logfile)) then
    call logfile%print_line('Fractional L2 error in unitarity of '// &
       & matrix_name//': '//error)
  endif
  
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' not as unitary as expected. &
       &Fractional L2 error in unitarity: '//error)
  endif
end subroutine

subroutine check_orthonormal_RealMatrix_character(input,matrix_name,logfile, &
   & warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  character(*),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  
  call check_orthonormal(input,str(matrix_name),logfile,warning_threshold)
end subroutine

subroutine check_orthonormal_ComplexMatrix_String(input,matrix_name,logfile, &
   & warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  type(String),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  if (size(input,1)<size(input,2)) then
    error = sum_squares( input*hermitian(input) &
                     & - cmplxmat(make_identity_matrix(size(input,1))))
  else
    error = sum_squares( hermitian(input)*input &
                     & - cmplxmat(make_identity_matrix(size(input,2))))
  endif
  error = sqrt(error)
  
  if (present(logfile)) then
    call logfile%print_line('Fractional L2 error in unitarity of '// &
       & matrix_name//': '//error)
  endif
  
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' not as unitary as expected. &
       &Fractional L2 error in unitarity: '//error)
  endif
end subroutine

subroutine check_orthonormal_ComplexMatrix_character(input,matrix_name, &
   & logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  character(*),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  
  call check_orthonormal(input,str(matrix_name),logfile,warning_threshold)
end subroutine

subroutine check_real_String(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  type(String),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  error = sqrt(sum_squares(aimag(input)))
  
  if (present(logfile)) then
    call logfile%print_line('L2 error in imaginary components of '// &
       & matrix_name//': '//error)
  endif
  
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' not as real as expected. &
       &L2 error in imaginary components: '//error)
  endif
end subroutine

subroutine check_real_character(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(ComplexMatrix), intent(in)              :: input
  character(*),        intent(in)              :: matrix_name
  type(OFile),         intent(inout), optional :: logfile
  real(dp),            intent(in),    optional :: warning_threshold
  
  call check_real(input,str(matrix_name),logfile,warning_threshold)
end subroutine

subroutine check_int_String(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  type(String),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  
  real(dp) :: threshold
  real(dp) :: error
  
  if (present(warning_threshold)) then
    threshold = warning_threshold
  else
    threshold = 1e-6_dp
  endif
  
  error = sqrt(sum_squares(input-nint(input)))
  
  if (present(logfile)) then
    call logfile%print_line('L2 error in difference from integer of &
       &components of '//matrix_name//': '//error)
  endif
  
  if (error>threshold) then
    call print_line(WARNING//': '//matrix_name//' not as close to integers as &
       &expected. L2 error in components: '//error)
  endif
end subroutine

subroutine check_int_character(input,matrix_name,logfile,warning_threshold)
  implicit none
  
  type(RealMatrix), intent(in)              :: input
  character(*),     intent(in)              :: matrix_name
  type(OFile),      intent(inout), optional :: logfile
  real(dp),         intent(in),    optional :: warning_threshold
  
  call check_int(input,str(matrix_name),logfile,warning_threshold)
end subroutine
end module

! ======================================================================
! Routines for finding the eigenvalues and eigenvectors of Hermitian,
!    symmetric and general matrices.
! ======================================================================
module caesar_hermitian_eigenstuff_module
  use caesar_foundations_module
  use caesar_io_module
  
  use caesar_lapack_wrapper_module
  use caesar_linear_algebra_module
  use caesar_phase_module
  use caesar_qr_decomposition_module
  implicit none
  
  private
  
  ! Classes.
  public :: SymmetricEigenstuff
  public :: HermitianEigenstuff
  public :: ComplexEigenstuff
  
  ! Diagonalisation routines.
  public :: diagonalise_symmetric
  public :: diagonalise_hermitian
  public :: diagonalise_complex
  
  ! --------------------------------------------------
  ! The eigen(values and vectors) of a matrix.
  ! --------------------------------------------------
  type SymmetricEigenstuff
    real(dp)              :: eval
    real(dp), allocatable :: evec(:)
  end type
  
  type HermitianEigenstuff
    real(dp)                 :: eval
    complex(dp), allocatable :: evec(:)
  end type
  
  type ComplexEigenstuff
    complex(dp)              :: eval
    complex(dp), allocatable :: left_evec(:)
    complex(dp), allocatable :: right_evec(:)
  end type
  
  ! --------------------------------------------------
  ! Interfaces.
  ! --------------------------------------------------
  interface diagonalise_symmetric
    module procedure diagonalise_symmetric_reals
    module procedure diagonalise_symmetric_RealMatrix
  end interface
  
  interface diagonalise_hermitian
    module procedure diagonalise_hermitian_complexes
    module procedure diagonalise_hermitian_ComplexMatrix
  end interface
  
  interface diagonalise_complex
    module procedure diagonalise_complex_complexes
    module procedure diagonalise_complex_ComplexMatrix
  end interface
contains

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a real, symmetric matrix.
! --------------------------------------------------
function diagonalise_symmetric_reals(input,basis) result(output)
  implicit none
  
  real(dp),         intent(in)           :: input(:,:)
  type(RealVector), intent(in), optional :: basis(:)
  type(SymmetricEigenstuff), allocatable :: output(:)
  
  ! dsyev variables.
  integer               :: n       ! The size of the input.
  real(dp), allocatable :: a(:,:)  ! The input. Also eigenvectors.
  real(dp), allocatable :: w(:)    ! The eigenvalues.
  real(dp), allocatable :: work(:) ! Workspace.
  integer               :: lwork   ! Workspace.
  integer               :: info    ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  if (present(basis)) then
    if (size(basis)==0) then
      allocate(output(0), stat=ialloc); call err(ialloc)
      return
    endif
    
    do i=1,size(basis)
      if (size(basis(i))/=size(input,1)) then
        call print_line(CODE_ERROR//': Trying to find the eigenavlues of a &
           &matrix in an incompatible basis.')
      endif
    enddo
  endif
  
  if (.not. present(basis)) then
    a = input
  else
    a = row_matrix(basis)
    a = dble(mat(a) * mat(input) * transpose(mat(a)))
  endif
  
  n = size(a,1)
  allocate( w(n),        &
          & work(3*n-1), &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork
  call dsyev( jobz  = 'V',  &
            & uplo  = 'U',  &
            & n     = n,    &
            & a     = a,    &
            & lda   = n,    &
            & w     = w,    &
            & work  = work, &
            & lwork = -1,   &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: dsyev error code: '//info)
    call err()
  endif
  lwork = nint(work(1))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call dsyev( jobz  = 'V',   &
            & uplo  = 'U',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & work  = work,  &
            & lwork = lwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: dsyev error code: '//info)
    call print_line('matrix:')
    call print_lines(mat(input))
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    if (.not. present(basis)) then
      output(i)%evec = a(:,i)
    else
      output(i)%evec = dble(sum(a(:,i)*basis))
    endif
  enddo
end function

function diagonalise_symmetric_RealMatrix(input,basis) result(output)
  implicit none
  
  type(RealMatrix), intent(in)           :: input
  type(RealVector), intent(in), optional :: basis(:)
  type(SymmetricEigenstuff), allocatable :: output(:)
  
  output = diagonalise_symmetric(dble(input), basis)
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of a complex, Hermitian matrix.
! --------------------------------------------------
function diagonalise_hermitian_complexes(input,basis) result(output)
  implicit none
  
  complex(dp),         intent(in)           :: input(:,:)
  type(ComplexVector), intent(in), optional :: basis(:)
  type(HermitianEigenstuff), allocatable    :: output(:)
  
  ! zheev variables.
  integer                  :: n        ! The size of the input.
  complex(dp), allocatable :: a(:,:)   ! Input matrix. Also eigenvectors.
  real(dp),    allocatable :: w(:)     ! Eigenvalues.
  complex(dp), allocatable :: work(:)  ! Workspace.
  integer                  :: lwork    ! Workspace.
  real(dp),    allocatable :: rwork(:) ! Workspace.
  integer                  :: info     ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  if (present(basis)) then
    if (size(basis)==0) then
      allocate(output(0), stat=ialloc); call err(ialloc)
      return
    endif
    
    do i=1,size(basis)
      if (size(basis(i))/=size(input,1)) then
        call print_line(CODE_ERROR//': Trying to find the eigenavlues of a &
           &matrix in an incompatible basis.')
      endif
    enddo
  endif
  
  if (.not. present(basis)) then
    a = input
  else
    a = row_matrix(basis)
    a = cmplx(mat(a) * mat(input) * hermitian(mat(a)))
  endif
  
  n = size(a,1)
  allocate( w(n),         &
          & rwork(3*n-2), &
          & work(3*n-1),  &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork
  call zheev( jobz  = 'V',   &
            & uplo  = 'U',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & work  = work,  &
            & lwork = -1,    &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zheev error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call zheev( jobz  = 'V',   &
            & uplo  = 'U',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & work  = work,  &
            & lwork = lwork, &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zheev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    
    if (.not. present(basis)) then
      output(i)%evec = a(:,i)
    else
      output(i)%evec = cmplx(sum(a(:,i)*basis))
    endif
  enddo
end function

function diagonalise_hermitian_ComplexMatrix(input,basis) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)           :: input
  type(ComplexVector), intent(in), optional :: basis(:)
  type(HermitianEigenstuff), allocatable    :: output(:)
  
  output = diagonalise_hermitian(cmplx(input), basis)
end function

! --------------------------------------------------
! Calculates the eigenvalues and eigenvectors of
!    a general square complex matrix.
! --------------------------------------------------
function diagonalise_complex_complexes(input,basis) result(output)
  implicit none
  
  complex(dp),         intent(in)           :: input(:,:)
  type(ComplexVector), intent(in), optional :: basis(:)
  type(ComplexEigenstuff), allocatable      :: output(:)
  
  ! zgeev variables.
  integer                  :: n        ! The size of the input.
  complex(dp), allocatable :: a(:,:)   ! The input.
  complex(dp), allocatable :: w(:)     ! Eigenvalues.
  complex(dp), allocatable :: vl(:,:)  ! Left eigenvectors.
  complex(dp), allocatable :: vr(:,:)  ! Right eigenvectors.
  complex(dp), allocatable :: work(:)  ! Workspace.
  integer                  :: lwork    ! Workspace.
  real(dp),    allocatable :: rwork(:) ! Workspace.
  integer                  :: info     ! Error code. 0 on success.
  
  ! Temporary variables.
  integer :: i,ialloc
  
  if (size(input,1)==0) then
    allocate(output(0), stat=ialloc); call err(ialloc)
    return
  endif
  
  if (size(input,1)/=size(input,2)) then
    call print_line(CODE_ERROR//': Trying to find the eigenvalues of a &
       &non-square matrix.')
    call err()
  endif
  
  if (present(basis)) then
    if (size(basis)==0) then
      allocate(output(0), stat=ialloc); call err(ialloc)
      return
    endif
    
    do i=1,size(basis)
      if (size(basis(i))/=size(input,1)) then
        call print_line(CODE_ERROR//': Trying to find the eigenavlues of a &
           &matrix in an incompatible basis.')
      endif
    enddo
  endif
  
  if (.not. present(basis)) then
    a = input
  else
    a = row_matrix(basis)
    a = cmplx(mat(a) * mat(input) * transpose(mat(a)))
  endif
  
  n = size(a,1)
  allocate( w(n),       &
          & vl(n,n),    &
          & vr(n,n),    &
          & rwork(2*n), &
          & work(2*n),  &
          & stat=ialloc); call err(ialloc)
  
  ! calculate optimal lwork.
  call zgeev( jobvl = 'V',   &
            & jobvr = 'V',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & vl    = vl,    &
            & ldvl  = n,     &
            & vr    = vr,    &
            & ldvr  = n,     &
            & work  = work,  &
            & lwork = -1,    &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zgeev error code: '//info)
    call err()
  endif
  lwork = nint(real(work(1)))
  deallocate(work, stat=ialloc); call err(ialloc)
  allocate(work(lwork), stat=ialloc); call err(ialloc)
  
  ! calculate eigenstuff
  call zgeev( jobvl = 'V',   &
            & jobvr = 'V',   &
            & n     = n,     &
            & a     = a,     &
            & lda   = n,     &
            & w     = w,     &
            & vl    = vl,    &
            & ldvl  = n,     &
            & vr    = vr,    &
            & ldvr  = n,     &
            & work  = work,  &
            & lwork = lwork, &
            & rwork = rwork, &
            & info  = info)
  if (info /= 0) then
    call print_line(ERROR//' in diagonalisation: zheev error code: '//info)
    call err()
  endif
  
  allocate(output(n), stat=ialloc); call err(ialloc)
  do i=1,n
    output(i)%eval = w(i)
    if (.not. present(basis)) then
      output(i)%left_evec = vl(:,i)
      output(i)%right_evec = vr(:,i)
    else
      output(i)%left_evec = cmplx(sum(vl(:,i)*basis))
      output(i)%right_evec = cmplx(sum(vr(:,i)*conjg(basis)))
    endif
  enddo
end function

function diagonalise_complex_ComplexMatrix(input,basis) result(output)
  implicit none
  
  type(ComplexMatrix), intent(in)           :: input
  type(ComplexVector), intent(in), optional :: basis(:)
  type(ComplexEigenstuff), allocatable      :: output(:)
  
  output = diagonalise_complex(cmplx(input), basis)
end function
end module

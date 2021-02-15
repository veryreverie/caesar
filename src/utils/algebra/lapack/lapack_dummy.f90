! ======================================================================
! A dummy module for when BLAS/LAPACK is not being linked.
! ======================================================================
module caesar_lapack_wrapper_module
  use caesar_foundations_module
  use caesar_io_module
  implicit none
  
  private
  
  ! Vector copying.
  public :: dcopy
  public :: zcopy
  
  ! Dot product.
  public :: ddot
    
  ! Multiplication between vectors and scalars.
  public :: dscal
  public :: zscal
  
  ! Complex norm.
  public :: dznrm2
  
  ! Linear least-squares.
  public :: dgels
    
  ! LU factorisation.
  public :: dgetrf
    
  ! Matrix inversion.
  public :: dgetri
  
  ! Diagonalisation.
  public :: dsyev
  public :: zheev
  public :: zgeev
  
  ! QR decomposition.
  public :: dgeqrf
  public :: zgeqrf
  public :: dorgqr
  public :: zungqr
  
  public :: LAPACK_LINKED
  
  logical, parameter :: LAPACK_LINKED = .false.
contains
! Copies a real vector. Equivalent to dy = dx.
subroutine dcopy(n,dx,incx,dy,incy)
  implicit none
  
  integer,  intent(in)  :: n     ! Length of vectors
  real(dp), intent(in)  :: dx(*) ! Input vector
  integer,  intent(in)  :: incx  ! Increment along dx
  real(dp), intent(out) :: dy(*) ! Output vector
  integer,  intent(in)  :: incy  ! Increment along dy
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Copies complex vector. Equivalent to zy = zx.
subroutine zcopy(n,zx,incx,zy,incy)
  implicit none
  
  integer,     intent(in)  :: n     ! Length of vectors
  complex(dp), intent(in)  :: zx(*) ! Input vector
  integer,     intent(in)  :: incx  ! Increment along zx
  complex(dp), intent(out) :: zy(*) ! Output vector
  integer,     intent(in)  :: incy  ! Increment along zy
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Real dot product. Returns dx.dy.
function ddot(n,dx,incx,dy,incy) result(output)
  implicit none
  
  integer,  intent(in) :: n     ! Length of vectors
  real(dp), intent(in) :: dx(*) ! First vector
  integer,  intent(in) :: incx  ! Increment along dx
  real(dp), intent(in) :: dy(*) ! Second vector
  integer,  intent(in) :: incy  ! Increment along dy
  real(dp)             :: output
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end function

! Multiplies real vector by real scalar. Equivalent to dx *= da.
subroutine dscal(n,da,dx,incx)
  implicit none
  
  integer,  intent(in)    :: n     ! Length of vector
  real(dp), intent(in)    :: da    ! Scalar
  real(dp), intent(inout) :: dx(*) ! Vector
  integer,  intent(in)    :: incx  ! Increment along dx
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Multiplies complex vector by complex scalar. Equivalent to zx *= za.
subroutine zscal(n,za,zx,incx)
  implicit none
  
  integer,     intent(in)    :: n     ! Length of vector
  complex(dp), intent(in)    :: za    ! Scalar
  complex(dp), intent(inout) :: zx(*) ! Vector
  integer,     intent(in)    :: incx  ! Increment along zx
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Complex norm. Returns sqrt(x.x).
function dznrm2(n,x,incx) result(output)
  implicit none
  
  integer,     intent(in) :: n      ! Length of vector
  complex(dp), intent(in) :: x(*)   ! Vector
  integer,     intent(in) :: incx   ! Increment along x
  real(dp)                :: output ! Result
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end function

! Minimises the least-squares fit l=(a.x-b)**2.
subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
  implicit none
  
  character(1), intent(in)    :: trans    ! n/t: if t, a is transposed.
  integer,      intent(in)    :: m        ! size(a,1)=size(b,1).
  integer,      intent(in)    :: n        ! size(a,2)=size(x,1).
  integer,      intent(in)    :: nrhs     ! size(b,2)=size(x,2).
  integer,      intent(in)    :: lda      ! The leading dimension of a.
  real(dp),     intent(inout) :: a(lda,*) ! The matrix a.
  integer,      intent(in)    :: ldb      ! The leading dimension of b.
  real(dp),     intent(inout) :: b(ldb,*) ! The matrix b.
  real(dp),     intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,      intent(in)    :: lwork    ! size(work).
  integer,      intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! In-place LU factorisation. Required for dgetri. LU factorises a.
subroutine dgetrf(m,n,a,lda,ipiv,info)
  implicit none
  
  integer,  intent(in)    :: m        ! size(a,1).
  integer,  intent(in)    :: n        ! size(a,2).
  integer,  intent(in)    :: lda      ! The leading dimension of a.
  real(dp), intent(inout) :: a(lda,*) ! The matrix a.
  integer,  intent(out)   :: ipiv(*)  ! Pivot indices. size=min(n,m).
  integer,  intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! In-place matrix inversion. Inverts a.
subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
  implicit none
  
  integer,  intent(in)    :: n        ! size(a,1)=size(a,2).
  integer,  intent(in)    :: lda      ! The leading dimension of a.
  real(dp), intent(inout) :: a(lda,*) ! The matrix to be inverted.
  integer,  intent(in)    :: ipiv(*)  ! Pivot indices. size=n.
  real(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,  intent(in)    :: lwork    ! size(work).
  integer,  intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Finds the eigenvalues of a real symmetric matrix.
subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
  implicit none
  
  character(1), intent(in)    :: jobz     ! N/V: if v, calculate eigenvecs.
  character(1), intent(in)    :: uplo     ! U/L: upper/lower triangle.
  integer,      intent(in)    :: n        ! The order of a.
  integer,      intent(in)    :: lda      ! The dimension of a.
  real(dp),     intent(inout) :: a(lda,*) ! Symmetric matrix.
  real(dp),     intent(out)   :: w(*)     ! Eigenvalues of a.
  real(dp),     intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,      intent(in)    :: lwork    ! The length of work.
  integer,      intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Finds the eigenvalues of a complex hermitian matrix.
subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
  implicit none
  
  character(1), intent(in)    :: jobz     ! N/V: if v, calculate eigenvecs.
  character(1), intent(in)    :: uplo     ! U/L: upper/lower triangle.
  integer,      intent(in)    :: n        ! The order of a.
  integer,      intent(in)    :: lda      ! The dimension of a.
  complex(dp),  intent(inout) :: a(lda,*) ! Hermitian matrix.
  real(dp),     intent(out)   :: w(*)     ! Eigenvalues of a.
  complex(dp),  intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,      intent(in)    :: lwork    ! The length of work.
  real(dp),     intent(out)   :: rwork(*) ! Working array.
  integer,      intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Find the eigenvalues of a general complex matrix.
subroutine zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork, &
   &info)
  implicit none
  
  character(1), intent(in)    :: jobvl      ! N/V: if v, calc. left evecs.
  character(1), intent(in)    :: jobvr      ! N/V: if v, calc. right evecs.
  integer,      intent(in)    :: n          ! The order of a.
  integer,      intent(in)    :: lda        ! The dimension of a.
  complex(dp),  intent(inout) :: a(lda,*)   ! Complex matrix.
  complex(dp),  intent(out)   :: w(*)       ! Eigenvalues of a.
  integer,      intent(in)    :: ldvl       ! The dimension of vl.
  complex(dp),  intent(out)   :: vl(ldvl,*) ! Left eigenvectors of a.
  integer,      intent(in)    :: ldvr       ! The dimension of vr.
  complex(dp),  intent(out)   :: vr(ldvr,*) ! Right eigenvectors of a.
  complex(dp),  intent(out)   :: work(*)    ! work(1) = optimal lwork.
  integer,      intent(in)    :: lwork      ! The length of work.
  real(dp),     intent(out)   :: rwork(*)   ! Working array.
  integer,      intent(out)   :: info       ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Find the QR decomposition of a real matrix.
! N.B. the output is given in terms of reflectors.
subroutine dgeqrf(m,n,a,lda,tau,work,lwork,info)
  implicit none
  
  integer,  intent(in)    :: m        ! size(a,1).
  integer,  intent(in)    :: n        ! size(a,2).
  integer,  intent(in)    :: lda      ! The leading dimension of a.
  real(dp), intent(inout) :: a(lda,*) ! Matrix a. Also R and reflectors.
  real(dp), intent(out)   :: tau(*)   ! Factors of reflectors.
  real(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,  intent(in)    :: lwork    ! The length of work.
  integer,  intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Find the QR decomposition of a complex matrix.
! N.B. the output is given in terms of reflectors.
subroutine zgeqrf(m,n,a,lda,tau,work,lwork,info)
  implicit none
  
  integer,     intent(in)    :: m        ! size(a,1).
  integer,     intent(in)    :: n        ! size(a,2).
  integer,     intent(in)    :: lda      ! The leading dimension of a.
  complex(dp), intent(inout) :: a(lda,*) ! Matrix a. Also R and reflectors.
  complex(dp), intent(out)   :: tau(*)   ! Factors of reflectors.
  complex(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,     intent(in)    :: lwork    ! The length of work.
  integer,     intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Convert reflectors into the Q matrix.
subroutine dorgqr(m,n,k,a,lda,tau,work,lwork,info)
  implicit none
  
  integer,  intent(in)    :: m        ! size(a,1).
  integer,  intent(in)    :: n        ! size(a,2).
  integer,  intent(in)    :: k        ! No. reflectors.
  integer,  intent(in)    :: lda      ! The leading dimension of a.
  real(dp), intent(inout) :: a(lda,*) ! The matrix a. Also the output.
  real(dp), intent(in)    :: tau(*)   ! Factors of reflectors.
  real(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,  intent(in)    :: lwork    ! The length of work.
  integer,  intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine

! Convert reflectors into the Q matrix.
subroutine zungqr(m,n,k,a,lda,tau,work,lwork,info)
  implicit none
  
  integer,     intent(in)    :: m        ! size(a,1).
  integer,     intent(in)    :: n        ! size(a,2).
  integer,     intent(in)    :: k        ! No. reflectors.
  integer,     intent(in)    :: lda      ! The leading dimension of a.
  complex(dp), intent(inout) :: a(lda,*) ! The matrix a. Also the output.
  complex(dp), intent(in)    :: tau(*)   ! Factors of reflectors.
  complex(dp), intent(out)   :: work(*)  ! work(1) = optimal lwork.
  integer,     intent(in)    :: lwork    ! The length of work.
  integer,     intent(out)   :: info     ! 0 on success.
  
  call print_line(ERROR//': Cannot perform linear algebra because Caesar has &
     &not been linked against BLAS/LAPACK. Please use the CMake flag &
     &-DLINK_TO_LAPACK:LOGICAL=true to link against LAPACK.')
  call err()
end subroutine
end module

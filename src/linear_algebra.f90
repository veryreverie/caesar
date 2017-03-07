! ------------------------------------------------------------
! Assorted linear algebra / vector algebra subroutines.
! Includes interfaces for BLAS and LAPACK routines.
! ------------------------------------------------------------
module linear_algebra
  use constants, only : dp
  implicit none
  
  ! The eigen(values and vectors) of a matrix
  ! degeneracy(i) = 0 if eigenvector(i) is non-degenerate
  !               = j if eigenvector(i) is degenerate,
  ! where j is an arbitrary but unique id
  type RealEigenstuff
    real(dp), allocatable :: evals(:)
    real(dp), allocatable :: evecs(:,:)
    integer,  allocatable :: degeneracy(:)
  end type
  
  type ComplexEigenstuff
    real(dp),    allocatable :: evals(:)
    complex(dp), allocatable :: evecs(:,:)
    integer,     allocatable :: degeneracy(:)
  end type

  interface size
    module procedure size_RealEigenstuff
    module procedure size_ComplexEigenstuff
  end interface
  
  interface calculate_eigenstuff
    module procedure calculate_RealEigenstuff
    module procedure calculate_ComplexEigenstuff
  end interface
  
  ! ----------------------------------------
  ! BLAS / LAPACK interface
  ! ----------------------------------------
  interface
    
    ! Copies a real vector. Equivalent to DY = DX
    pure subroutine dcopy(N,DX,INCX,DY,INCY)
      use constants, only : dp
      implicit none
      
      integer,  intent(in)  :: N     ! length of vectors
      real(dp), intent(in)  :: DX(*) ! input vector
      integer,  intent(in)  :: INCX  ! increment along DX
      real(dp), intent(out) :: DY(*) ! output vector
      integer,  intent(in)  :: INCY  ! increment along DY
    end subroutine
    
    ! Copies complex vector. Equivalent to ZY = ZX
    pure subroutine zcopy(N,ZX,INCX,ZY,INCY)
      use constants, only : dp
      implicit none
      
      integer,     intent(in)  :: N     ! length of vectors
      complex(dp), intent(in)  :: ZX(*) ! input vector
      integer,     intent(in)  :: INCX  ! increment along ZX
      complex(dp), intent(out) :: ZY(*) ! output vector
      integer,     intent(in)  :: INCY  ! increment along ZY
    end subroutine
    
    ! Real dot product. Returns DX.DY
    pure real(dp) function ddot(N,DX,INCX,DY,INCY)
      use constants, only : dp
      implicit none
      
      integer,  intent(in) :: N     ! length of vectors
      real(dp), intent(in) :: DX(*) ! first vector
      integer,  intent(in) :: INCX  ! increment along DX
      real(dp), intent(in) :: DY(*) ! second vector
      integer,  intent(in) :: INCY  ! increment along DY
    end function
    
    ! Multiplies real vector by real scalar. Equivalent to DX *= DA
    pure subroutine dscal(N,DA,DX,INCX)
      use constants, only : dp
      implicit none
      
      integer,  intent(in)    :: N     ! length of vector
      real(dp), intent(in)    :: DA    ! scalar
      real(dp), intent(inout) :: DX(*) ! vector
      integer,  intent(in)    :: INCX  ! increment along DX
    end subroutine
    
    ! Multiplies complex vector by complex scalar. Equivalent to ZX *= ZA
    pure subroutine zscal(N,ZA,ZX,INCX)
      use constants, only : dp
      implicit none
      
      integer,     intent(in)    :: N     ! length of vector
      complex(dp), intent(in)    :: ZA    ! scalar
      complex(dp), intent(inout) :: ZX(*) ! vector
      integer,     intent(in)    :: INCX  ! increment along ZX
    end subroutine
    
    ! Complex norm. Returns sqrt(X.X)
    pure function dznrm2(N,X,INCX) result(output)
      use constants, only : dp
      implicit none
      
      integer,     intent(in) :: N      ! length of vector
      complex(dp), intent(in) :: X(*)   ! vector
      integer,     intent(in) :: INCX   ! increment along X
      real(dp)                :: output ! result
    end function
    
    ! Finds the eigenvalues of a hermitian matrix
    pure subroutine zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
      use constants, only : dp
      implicit none
      
      character(1), intent(in)    :: JOBZ     ! N/V: if V, calculate eigenvectors
      character(1), intent(in)    :: UPLO     ! U/L: store upper/lower triangle
      integer,      intent(in)    :: N        ! the order of A
      integer,      intent(in)    :: LDA      ! the dimension of A
      complex(dp),  intent(inout) :: A(LDA,*) ! Hermitian matrix
      real(dp),     intent(out)   :: W(*)     ! eigenvalues of A
      complex(dp),  intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK
      integer,      intent(in)    :: LWORK    ! the length of WORK
      real(dp),     intent(out)   :: RWORK(*) ! working array
      integer,      intent(out)   :: INFO     ! 0 on success
    end subroutine
    
    ! Finds the eigenvalues of a symmetric matrix
    pure subroutine dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
      use constants, only : dp
      implicit none
      
      character(1), intent(in)    :: JOBZ     ! N/V: if V, calculate eigenvectors
      character(1), intent(in)    :: UPLO     ! U/L: store upper/lower triangle
      integer,      intent(in)    :: N        ! the order of A
      integer,      intent(in)    :: LDA      ! the dimension of A
      real(dp),     intent(inout) :: A(LDA,*) ! symmetric matrix
      real(dp),     intent(out)   :: W(*)     ! eigenvalues of A
      real(dp),     intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK
      integer,      intent(in)    :: LWORK    ! the length of WORK
      integer,      intent(out)   :: INFO     ! 0 on success
    end subroutine
  end interface

  ! ----------------------------------------
  ! determinant interface
  ! ----------------------------------------
  interface determinant
    module procedure determinant_integer, determinant_real
  end interface

contains

! ----------------------------------------
! Returns the number of states of an Eigenstuff
! ----------------------------------------
pure function size_RealEigenstuff(estuff) result(output)
  implicit none
  
  type(RealEigenstuff), intent(in) :: estuff
  integer                          :: output
  
  output = size(estuff%evals)
end function

pure function size_ComplexEigenstuff(estuff) result(output)
  implicit none
  
  type(ComplexEigenstuff), intent(in) :: estuff
  integer                             :: output
  
  output = size(estuff%evals)
end function

! ----------------------------------------
! given a 3x3 matrix A, returns det(A)
! ----------------------------------------
function determinant_integer(A) result(determinant)
  implicit none
  
  integer, intent(in) :: A(3,3)
  integer             :: determinant
  
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
             &+ A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))&
             &+ A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end function

function determinant_real(A) result(determinant)
  implicit none
  
  real(dp), intent(in) :: A(3,3)
  real(dp)             :: determinant
  
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
             &+ A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))&
             &+ A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end function

! ----------------------------------------
! calculates the inverse, B, of matrix A
! A and B are real, 3x3 matrices
! ----------------------------------------
function invert(A) result(B)
  use string_module
  implicit none
  
  real(dp), intent(in)  :: A(3,3)
  real(dp)              :: B(3,3)
  
  real(dp) :: d      ! 1/det(A)
  real(dp) :: C(3,3) ! transpose(A)
  
  d = 1.0_dp/determinant(A)
  
  ! check for d=infinity or d=NaN
  if (dabs(d)>huge(0.0_dp) .or. d<d) then
    call print_line('Error in invert: singular matrix.')
    stop
  endif
  
  C = transpose(A)
  
  B(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))*d
  B(1,2) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))*d
  B(1,3) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))*d
  B(2,1) = (C(3,2)*C(1,3)-C(3,3)*C(1,2))*d
  B(2,2) = (C(3,3)*C(1,1)-C(3,1)*C(1,3))*d
  B(2,3) = (C(3,1)*C(1,2)-C(3,2)*C(1,1))*d
  B(3,1) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))*d
  B(3,2) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))*d
  B(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))*d
end function

! ----------------------------------------
! Calculates B=inverse(A)*|det(A)|
! A and B are 3x3 integer matrices
! ----------------------------------------
function invert_int(A) result(B)
  implicit none
  
  integer, intent(in)  :: A(3,3)
  integer              :: B(3,3)
  
  integer :: C(3,3) ! transpose(A)
  integer :: d      ! 1 if det(A)>=0, -1 otherwise
  
  d = sign(1,determinant(A))
  
  C = transpose(A)
  
  B(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))*d
  B(1,2) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))*d
  B(1,3) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))*d
  B(2,1) = (C(3,2)*C(1,3)-C(3,3)*C(1,2))*d
  B(2,2) = (C(3,3)*C(1,1)-C(3,1)*C(1,3))*d
  B(2,3) = (C(3,1)*C(1,2)-C(3,2)*C(1,1))*d
  B(3,1) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))*d
  B(3,2) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))*d
  B(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))*d
end function

! Calculates the eigenvalues and eigenvectors of a real, symmetric matrix
function calculate_RealEigenstuff(input) result(output)
  use string_module
  implicit none
  
  real(dp), intent(in) :: input(:,:)  ! a real, symmetric matrix
  type(RealEigenstuff) :: output      ! the eigenvalues and eigenstates of a
  
  ! working variables
  integer               :: n
  real(dp), allocatable :: work(:)
  integer               :: lwork
  integer               :: info
  
  n = size(input,1)
  allocate(output%evals(n))
  allocate(output%evecs(n,n))
  output%evecs = input
  
  ! calculate optimal lwork
  allocate(work(3*n-1))
  call dsyev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), -1, info)
  if (info /= 0) then
    call print_line("dsyev failed, info= "//str(info))
    stop
  endif
  lwork = nint(work(1))
  deallocate(work)
  allocate(work(lwork))
  
  ! calculate eigenstuff
  call dsyev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), lwork, info)
  if (info /= 0) then
    call print_line("dsyev failed, info= "//str(info))
    stop
  endif
end function

! Calculates the eigenvalues and eigenvectors of a complex, hermitian matrix
function calculate_ComplexEigenstuff(input) result(output)
  use string_module
  implicit none
  
  complex(dp), intent(in) :: input(:,:)  ! a complex, hermitian matrix
  type(ComplexEigenstuff) :: output      ! the eigenvalues and eigenstates of a
  
  ! working variables
  integer               :: n
  complex(dp), allocatable :: work(:)
  integer               :: lwork
  real(dp), allocatable :: rwork(:)
  integer               :: info
  
  n = size(input,1)
  allocate(output%evals(n))
  allocate(output%evecs(n,n))
  output%evecs = input
  
  allocate(rwork(3*n-2))
  
  ! calculate optimal lwork
  allocate(work(3*n-1))
  call zheev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), -1, rwork, info)
  if (info /= 0) then
    call print_line("dsyev failed, info= "//str(info))
    stop
  endif
  lwork = nint(real(work(1)))
  deallocate(work)
  allocate(work(lwork))
  
  ! calculate eigenstuff
  call zheev('V', 'U', n, output%evecs(1,1), n, output%evals, &
    & work(1), lwork, rwork, info)
  if (info /= 0) then
    call print_line("zheev failed, info= "//str(info))
    stop
  endif
end function
end module

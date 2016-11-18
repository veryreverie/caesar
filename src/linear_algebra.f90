! ------------------------------------------------------------
! Assorted linear algebra / vector algebra subroutines.
! Includes interfaces for BLAS and LAPACK routines.
! ------------------------------------------------------------
module linear_algebra
  use constants, only : dp
  implicit none

! ----------------------------------------
! BLAS / LAPACK interface
! ----------------------------------------
interface
  ! Real dot product. Returns DX.DY
  real(dp) function ddot(N,DX,INCX,DY,INCY)
    use constants, only : dp
    implicit none
    
    integer,  intent(in) :: N     ! length of vectors
    real(dp), intent(in) :: DX(*) ! first vector
    integer,  intent(in) :: INCX  ! increment along DX
    real(dp), intent(in) :: DY(*) ! second vector
    integer,  intent(in) :: INCY  ! increment along DY
  end function
  
  ! Multiplies real vector by real scalar. Equivalent to DX *= DA
  subroutine dscal(N,DA,DX,INCX)
    use constants, only : dp
    implicit none
    
    integer,  intent(in)    :: N     ! length of vector
    real(dp), intent(in)    :: DA    ! scalar
    real(dp), intent(inout) :: DX(*) ! vector
    integer,  intent(in)    :: INCX  ! increment along DX
  end subroutine
  
  ! Multiplies complex vector by complex scalar. Equivalent to ZX *= ZA
  subroutine zscal(N,ZA,ZX,INCX)
    use constants, only : dp
    implicit none
    
    integer,     intent(in)    :: N     ! length of vector
    complex(dp), intent(in)    :: ZA    ! scalar
    complex(dp), intent(inout) :: ZX(*) ! vector
    integer,     intent(in)    :: INCX  ! increment along ZX
  end subroutine
  
  ! Copies complex vector. Equivalent to ZY = ZX
  subroutine zcopy(N,ZX,INCX,ZY,INCY)
    use constants, only : dp
    implicit none
    
    integer,     intent(in)  :: N     ! length of vectors
    complex(dp), intent(in)  :: ZX(*) ! input vector
    integer,     intent(in)  :: INCX  ! increment along ZX
    complex(dp), intent(out) :: ZY(*) ! output vector
    integer,     intent(in)  :: INCY  ! increment along ZY
  end subroutine
  
  ! Complex norm. Returns sqrt(X.X)
  real(kind(1.d0)) function dznrm2(N,X,INCX)
    use constants, only : dp
    implicit none
    
    integer,     intent(in) :: N    ! length of vector
    complex(dp), intent(in) :: X(*) ! vector
    integer,     intent(in) :: INCX ! increment along X
  end function
  
  ! Finds the eigenvalues of a hermitian matrix
  subroutine zheev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
    use constants, only : dp
    implicit none
    
    character(1), intent(in)    :: JOBZ     ! N/V: if V, calculate eigenvectors
    character(1), intent(in)    :: UPLO     ! U/L: store upper/lower triangle
    integer,      intent(in)    :: N        ! the order of A
    complex(dp),  intent(inout) :: A(LDA,*) ! Hermitian matrix
    integer,      intent(in)    :: LDA      ! the dimension of A
    real(dp),     intent(out)   :: W(*)     ! eigenvalues of A
    complex(dp),  intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK
    integer,      intent(in)    :: LWORK    ! the length of WORK
    real(dp),     intent(out)   :: RWORK(*) ! working array
    integer,      intent(out)   :: INFO     ! 0 on success
  end subroutine
  
  ! Finds the eigenvalues of a symmetric matrix
  subroutine dsyev(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
    use constants, only : dp
    implicit none
    
    character(1), intent(in)    :: JOBZ     ! N/V: if V, calculate eigenvectors
    character(1), intent(in)    :: UPLO     ! U/L: store upper/lower triangle
    integer,      intent(in)    :: N        ! the order of A
    real(dp),     intent(inout) :: A(LDA,*) ! symmetric matrix
    integer,      intent(in)    :: LDA      ! the dimension of A
    real(dp),     intent(out)   :: W(*)     ! eigenvalues of A
    real(dp),     intent(out)   :: WORK(*)  ! WORK(1) = optimal LWORK
    integer,      intent(in)    :: LWORK    ! the length of WORK
    integer,      intent(out)   :: INFO     ! 0 on success
  end subroutine
end interface

! ----------------------------------------
! determinant interface
! ----------------------------------------
interface determinant33
  module procedure determinant33_integer, determinant33_real
end interface

contains

! ----------------------------------------
! given a 3x3 matrix A, returns det(A)
! ----------------------------------------
function determinant33_integer(A) result(determinant)
  implicit none
  
  integer, intent(in) :: A(3,3)
  integer             :: determinant
  
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
             &+ A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
             &+ A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
end function

function determinant33_real(A) result(determinant)
  implicit none
  
  real(dp), intent(in) :: A(3,3)
  real(dp)             :: determinant
  
  determinant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
             &+ A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
             &+ A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
end function

! ----------------------------------------
! calculates the inverse, B, of matrix A
! A and B are real, 3x3 matrices
! ----------------------------------------
! TODO: |d|<epsilon would be a more stable check than d==0
subroutine inv_33(A,B)
  implicit none
  
  real(dp), intent(in)  :: A(3,3)
  real(dp), intent(out) :: B(3,3)
  real(dp)              :: d      ! det(A)
  
  d = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
   &+ A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
   &+ A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
  
  if (d==0.d0) then
    write(*,*) 'Error in inv_33: singular matrix.'
    stop
  endif
  
  d = 1.d0/d
  
  B(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))*d
  B(1,2) = (A(3,2)*A(1,3)-A(1,3)*A(3,2))*d
  B(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(3,2))*d
  B(2,1) = (A(3,1)*A(2,3)-A(2,3)*A(3,2))*d
  B(2,2) = (A(1,1)*A(3,3)-A(2,3)*A(3,2))*d
  B(2,3) = (A(2,1)*A(1,3)-A(2,3)*A(3,2))*d
  B(3,1) = (A(2,1)*A(3,2)-A(2,3)*A(3,2))*d
  B(3,2) = (A(3,1)*A(1,2)-A(2,3)*A(3,2))*d
  B(3,3) = (A(1,1)*A(2,2)-A(2,3)*A(3,2))*d
end subroutine

end module

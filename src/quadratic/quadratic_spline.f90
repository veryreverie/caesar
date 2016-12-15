!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Program file name: spline.f90                                          !
!                                                                         !
!  © Tao Pang 2006                                                        !
!                                                                         !
!  Last modified: August 1, 2012                                          !
!                                                                         !
!  (1) This F90 program is created for the book, "An Introduction to      !
!      Computational Physics, 2nd Edition," written by Tao Pang and       !
!      published by Cambridge University Press on January 19, 2006.       !
!                                                                         !
!  (2) No warranties, express or implied, are made for this program.      !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
module quadratic_spline_module
  implicit none
contains

subroutine quadratic_spline()
  use constants, only : dp
! An example of creating cubic-spline approximation of
! a discrete function fi=f(xi).
!
  IMPLICIT NONE
  INTEGER :: N, M, P
  INTEGER :: I, K
  REAL(dp) :: X, F, DX, H, ALPHA, BETA, GAMMA, ETA
  REAL(dp),ALLOCATABLE  :: XI(:), FI(:), P2(:)
  ! For reading input data
  INTEGER :: j

  ! Read in number of integration points
  OPEN(1,FILE='fit_input.dat')
  READ(1,*)N,M
  CLOSE(1)
  !N=N-1
 
  ALLOCATE(XI(N+1),FI(N+1),P2(N+1))

  OPEN(1,FILE='fit_energy.dat')
  DO j=1,N+1
    READ(1,*)XI(j),FI(j)
  !  WRITE(*,*)XI(j),FI(j)
  ENDDO ! N


  CALL CUBIC_SPLINE(N, XI, FI, P2)

!
! Find the approximation of the function
!
  
  
  OPEN(1,FILE='indep_pot.dat')
  write(1,*)XI(1),FI(1)
  H = (XI(N+1)-XI(1))/M 
  X = XI(1)
  
  P=M-1
  
  DO I = 1, P
      X = X + H
!
! Find the interval that x resides
    K = 1
    DX = X-XI(1)
    DO WHILE (DX .GE. 0)
      K = K + 1
      DX = X-XI(K)
    END DO
    K = K - 1
!
! Find the value of function f(x)
    DX = XI(K+1) - XI(K)
    ALPHA = P2(K+1)/(6*DX)
    BETA = -P2(K)/(6*DX)
    GAMMA = FI(K+1)/DX - DX*P2(K+1)/6
    ETA = DX*P2(K)/6 - FI(K)/DX
    F = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
       +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
       +GAMMA*(X-XI(K))+ETA*(X-XI(K+1))
    WRITE (1, *) X, F
  END DO
  CLOSE(1)

end subroutine


! as above, but takes arguments rather than reading files
! n.b. n is 'N+1' above.
! An example of creating cubic-spline approximation of
! a discrete function fi=f(xi).
pure function quadratic_spline2(m,xifi) result(output)
  use constants, only : dp
  implicit none
  
  integer,               intent(in) :: m           ! no. points to spline to
  real(dp), allocatable, intent(in) :: xifi(:,:)   ! {(x,f)} at n points
  real(dp), allocatable             :: output(:,:) ! {(x,f)} at m points
  
  integer :: n ! = size(xi) = size(xifi,2)
  real(dp), allocatable :: xi(:)
  real(dp), allocatable :: fi(:)
  
  integer :: p
  integer :: i, k
  real(dp) :: x, f, dx, h, alpha, beta, gamma, eta
  real(dp),allocatable  :: p2(:)

  n = size(xifi,2)
  
  allocate(xi(n))
  allocate(fi(n))
  allocate(p2(n))
  allocate(output(2,m))
  
  xi = xifi(1,:)
  fi = xifi(2,:)

  call cubic_spline(n-1, xi, fi, p2)

  ! find the approximation of the function
  
  output(1,1) = xi(1)
  output(2,1) = fi(1)
  
  h = (xi(n)-xi(1))/m 
  x = xi(1)
  
  p=m-1
  
  do i = 1, p
    x = x + h
    ! find the interval that x resides
    k = 1
    dx = x-xi(1)
    do while (dx .ge. 0)
      k = k + 1
      dx = x-xi(k)
    end do
    k = k - 1
  
    ! find the value of function f(x)
    dx = xi(k+1) - xi(k)
    alpha = p2(k+1)/(6*dx)
    beta = -p2(k)/(6*dx)
    gamma = fi(k+1)/dx - dx*p2(k+1)/6
    eta = dx*p2(k)/6 - fi(k)/dx
    f = alpha*(x-xi(k))*(x-xi(k))*(x-xi(k)) &
       +beta*(x-xi(k+1))*(x-xi(k+1))*(x-xi(k+1)) &
       +gamma*(x-xi(k))+eta*(x-xi(k+1))
    output(1,i+1) = x
    output(2,i+1) = f
  end do

end function


pure SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
  use constants, only : dp
  implicit none
  
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  INTEGER :: I
  INTEGER, INTENT (IN) :: N
  REAL(dp), INTENT (IN), DIMENSION (N+1):: XI, FI
  REAL(dp), INTENT (OUT), DIMENSION (N+1):: P2
  REAL(dp), DIMENSION (N):: G, H
  REAL(dp), DIMENSION (N-1):: D, B, C
!
! Assign the intervals and function differences
!
  DO I = 1, N
    H(I) = XI(I+1) - XI(I)
    G(I) = FI(I+1) - FI(I)
  END DO
!
! Evaluate the coefficient matrix elements
  DO I = 1, N-1
    D(I) = 2*(H(I+1)+H(I))
    B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
    C(I) = H(I+1)
  END DO
!
! Obtain the second-order derivatives
!
  CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
  P2(1) = 0
  P2(N+1) = 0
  DO I = 2, N 
    P2(I) = G(I-1)
  END DO
END SUBROUTINE CUBIC_SPLINE
!
pure SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
  use constants, only : dp
  implicit none
  
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  REAL(dp), INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL(dp), INTENT (OUT), DIMENSION (L):: Z
  REAL(dp), DIMENSION (L):: Y, W
  REAL(dp), DIMENSION (L-1):: V, T
!
! Evaluate the elements in the LU decomposition
!
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO
END SUBROUTINE TRIDIAGONAL_LINEAR_EQ
end module

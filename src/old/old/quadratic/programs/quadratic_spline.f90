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
PROGRAM SPLINE

! An example of creating cubic-spline approximation of
! a discrete function fi=f(xi).
!
  IMPLICIT NONE
  INTEGER :: N, M, P
  INTEGER :: I, K
  REAL :: X, F, DX, H, ALPHA, BETA, GAMMA, ETA
  REAL,ALLOCATABLE  :: XI(:), FI(:), P2(:)
  ! For reading input data
  INTEGER :: ierr,j
  CHARACTER(200) :: char200
  CHARACTER(80) :: in_file,mode

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

END PROGRAM SPLINE
!

SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  INTEGER :: I
  INTEGER, INTENT (IN) :: N
  REAL, INTENT (IN), DIMENSION (N+1):: XI, FI
  REAL, INTENT (OUT), DIMENSION (N+1):: P2
  REAL, DIMENSION (N):: G, H
  REAL, DIMENSION (N-1):: D, B, C
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
SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  REAL, INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL, INTENT (OUT), DIMENSION (L):: Z
  REAL, DIMENSION (L):: Y, W
  REAL, DIMENSION (L-1):: V, T
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


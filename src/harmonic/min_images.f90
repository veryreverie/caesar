! Subroutines for the calculation of minimum-image distances.
MODULE min_images
  USE utils,ONLY : dp,errstop
  IMPLICIT NONE
  
  ! Maximum number of images
  integer, parameter :: maxim = 8
  
interface min_images_brute_force
  module procedure min_images_brute_force_a ! takes rec_vec
  module procedure min_images_brute_force_b ! does not take rec_vec
end interface
  
CONTAINS

! This subroutine computes the minimum image vector(s) b of
! vector a with respect to the lattice specified by the columns of 
! lat_vec.  rec_vec are the reciprocal lattice vectors (w/o 2pi).
! -b is the vector from a to its closest lattice point.  nim is the number
! of image vectors.
subroutine min_images_brute_force_a(a,lat_vec,rec_vec,b,nim)
  implicit none
  
  REAL(dp),INTENT(in) :: a(3),lat_vec(3,3),rec_vec(3,3)
  REAL(dp),INTENT(out) :: b(3,maxim)
  INTEGER,INTENT(out) :: nim
  
  REAL(dp) :: Delta1(3),Delta2(3),Delta3(3),mag_b_sq,dist2,tol_L2
  INTEGER :: n(3),i,j,k
  ! Number of "shells" of lattice points to check.  Only used in setup, so
  ! may as well overkill.
  INTEGER,PARAMETER :: check_shell=3
  REAL(dp),PARAMETER :: tol=1.d-9
  tol_L2=tol*DOT_PRODUCT(lat_vec(1:3,1),lat_vec(1:3,1))
  n(1)=FLOOR(DOT_PRODUCT(a(1:3),rec_vec(1:3,1)))
  n(2)=FLOOR(DOT_PRODUCT(a(1:3),rec_vec(1:3,2)))
  n(3)=FLOOR(DOT_PRODUCT(a(1:3),rec_vec(1:3,3)))
  mag_b_sq=-1.d0
  nim=-1
  DO i=n(1)-check_shell,n(1)+check_shell+1
    Delta1=a-DBLE(i)*lat_vec(1:3,1)
    DO j=n(2)-check_shell,n(2)+check_shell+1
      Delta2=Delta1-DBLE(j)*lat_vec(1:3,2)
      DO k=n(3)-check_shell,n(3)+check_shell+1
        Delta3=Delta2-DBLE(k)*lat_vec(1:3,3)
        dist2=DOT_PRODUCT(Delta3,Delta3)
        IF(ABS(dist2-mag_b_sq)<=tol_L2)THEN
          nim=nim+1
          IF(nim>maxim)CALL errstop('MIN_IMAGES_BRUTE_FORCE', &
            &'Need to increase maxim parameter.')
          b(1:3,nim)=Delta3(1:3)
        ELSEIF(dist2<mag_b_sq.OR.nim==-1)THEN
          mag_b_sq=dist2
          nim=1
          b(1:3,1)=Delta3(1:3)
        ENDIF
      ENDDO ! k
    ENDDO ! j
  ENDDO ! i
  IF(nim<=0)CALL errstop('MIN_IMAGES_BRUTE_FORCE','Bug.')
end subroutine


! Compute the minimum image vector(s) b of vector a with respect to the 
! lattice specified by the columns of lat_vec. rec_vec are the reciprocal 
! lattice vectors (w/o 2pi). -b is the vector from a to its closest lattice 
! point.  nim is the number of image vectors.
SUBROUTINE min_images_brute_force_b(a,lat_vec,b,nim)
  use linear_algebra, only : invert
  implicit none
  
  REAL(dp),INTENT(in) :: a(3),lat_vec(3,3)
  REAL(dp),INTENT(out) :: b(3,8)
  INTEGER,INTENT(out) :: nim
  
  REAL(dp) :: rec_vec(3,3),delta1(3),delta2(3),delta3(3),mag_b_sq,dist2,tol_L2
  INTEGER :: n(3),i,j,k
  INTEGER,PARAMETER :: check_shell=3
  REAL(dp),PARAMETER :: tol=1.d-8
  
  rec_vec = transpose(invert(lat_vec))
  
  tol_L2=tol*dot_product(lat_vec(1,1:3),lat_vec(1,1:3))
  n(1)=floor(dot_product(a(1:3),rec_vec(1,1:3)))
  n(2)=floor(dot_product(a(1:3),rec_vec(2,1:3)))
  n(3)=floor(dot_product(a(1:3),rec_vec(3,1:3)))
  
  mag_b_sq=-1.d0
  nim=-1
  
  do i=n(1)-check_shell,n(1)+check_shell+1
    delta1=a-dble(i)*lat_vec(1,1:3)
    do j=n(2)-check_shell,n(2)+check_shell+1
      delta2=delta1-dble(j)*lat_vec(2,1:3)
      do k=n(3)-check_shell,n(3)+check_shell+1
        delta3=delta2-dble(k)*lat_vec(3,1:3)
        dist2=dot_product(delta3,delta3)
        if(abs(dist2-mag_b_sq)<=tol_L2)then
          nim=nim+1
          if(nim>8)call errstop('MIN_IMAGES_BRUTE_FORCE','Need to increase &
            &maxim parameter.')
          b(1:3,nim)=delta3(1:3)
        elseif(dist2<mag_b_sq.or.nim==-1)then
          mag_b_sq=dist2
          nim=1
          b(1:3,1)=delta3(1:3)
        endif
      enddo ! k
    enddo ! j
  enddo ! i
  
  if(nim<=0)call errstop('MIN_IMAGES_BRUTE_FORCE','Bug.')
end subroutine
END MODULE min_images

MODULE min_images
  ! Subroutines for the calculation of minimum-image distances.
  USE utils,ONLY : dp,errstop
  IMPLICIT NONE
  ! Maximum possible number of images.
  INTEGER,PARAMETER :: maxim=8


CONTAINS


  SUBROUTINE min_images_brute_force(a,lat_vec,rec_vec,b,nim)
    ! This subroutine computes the minimum image vector(s) b of
    ! vector a with respect to the lattice specified by the columns of 
    ! lat_vec.  rec_vec are the reciprocal lattice vectors (w/o 2pi).
    ! -b is the vector from a to its closest lattice point.  nim is the number
    ! of image vectors.
    IMPLICIT NONE
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
  END SUBROUTINE min_images_brute_force


  LOGICAL FUNCTION is_lat_point(rvec,rec_vec)
    ! This function returns T if and only if rvec is a lattice vector.  rec_vec
    ! holds the reciprocal lattice vectors.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: rvec(3),rec_vec(3,3)
    REAL(dp) :: t
    REAL(dp),PARAMETER :: tol=1.d-3
    t=DOT_PRODUCT(rvec,rec_vec(1:3,1))
    IF(ABS(ANINT(t)-t)<tol)THEN
      t=DOT_PRODUCT(rvec,rec_vec(1:3,2))
      IF(ABS(ANINT(t)-t)<tol)THEN
        t=DOT_PRODUCT(rvec,rec_vec(1:3,3))
        is_lat_point=(ABS(ANINT(t)-t)<tol)
      ELSE
        is_lat_point=.FALSE.
      ENDIF ! 2nd component integer
    ELSE
      is_lat_point=.FALSE.
    ENDIF ! 1st component integer
  END FUNCTION is_lat_point


END MODULE min_images
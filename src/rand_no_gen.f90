MODULE rand_no_gen
  ! Pseudo-random number generator.
  USE utils,ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC ranx
  INTEGER :: iseed=-1 ! Seed.  Supply a negative integer.


CONTAINS


  REAL(dp) FUNCTION ranx()
    ! Random number generator, adapted from ran2 in Numerical Recipes.
    ! (Method of l'Ecuyer with Bays-Durham shuffle.)
    IMPLICIT NONE
    INTEGER,PARAMETER :: im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014, &
      &ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32, &
      &ndiv=1+imm1/ntab,ntabp8=ntab+8
    INTEGER :: j,k
    INTEGER,SAVE :: iseed2=123456789,iv(ntab)=0,iy=0
    REAL(dp),PARAMETER :: am=1.d0/im1,rnmx=1.d0-EPSILON(1.d0)
    IF(iseed<=0)THEN
      iseed=MAX(-iseed,1)
      iseed2=iseed
      DO j=ntabp8,1,-1
        k=iseed/iq1
        iseed=ia1*(iseed-k*iq1)-k*ir1
        IF(iseed<0)iseed=iseed+im1
        IF(j<=ntab)iv(j)=iseed
      ENDDO ! j
      iy=iv(1)
    ENDIF ! iseed<=0
    k=iseed/iq1
    iseed=ia1*(iseed-k*iq1)-k*ir1
    IF(iseed<0)iseed=iseed+im1
    k=iseed2/iq2
    iseed2=ia2*(iseed2-k*iq2)-k*ir2
    IF(iseed2<0)iseed2=iseed2+im2
    j=1+iy/ndiv
    iy=iv(j)-iseed2
    iv(j)=iseed
    IF(iy<1)iy=iy+imm1
    ranx=MIN(am*iy,rnmx)
  END FUNCTION ranx


END MODULE rand_no_gen

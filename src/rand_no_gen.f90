! Pseudo-random number generator.
module rand_no_gen
  use utils, only : dp
  implicit none
  
  private
  public ranx
  integer :: iseed=-1 ! seed.  supply a negative integer.

contains

! ----------------------------------------------------------------------
! Random number generator, adapted from ran2 in Numerical Recipes.
! (Method of l'Ecuyer with Bays-Durham shuffle.)
! ----------------------------------------------------------------------
real(dp) function ranx()
  implicit none
  integer,parameter :: im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014, &
    &ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab=32, &
    &ndiv=1+imm1/ntab,ntabp8=ntab+8
  integer :: j,k
  integer,save :: iseed2=123456789,iv(ntab)=0,iy=0
  real(dp),parameter :: am=1.d0/im1,rnmx=1.d0-epsilon(1.d0)
  if(iseed<=0)then
    iseed=max(-iseed,1)
    iseed2=iseed
    do j=ntabp8,1,-1
      k=iseed/iq1
      iseed=ia1*(iseed-k*iq1)-k*ir1
      if(iseed<0)iseed=iseed+im1
      if(j<=ntab)iv(j)=iseed
    enddo ! j
    iy=iv(1)
  endif ! iseed<=0
  k=iseed/iq1
  iseed=ia1*(iseed-k*iq1)-k*ir1
  if(iseed<0)iseed=iseed+im1
  k=iseed2/iq2
  iseed2=ia2*(iseed2-k*iq2)-k*ir2
  if(iseed2<0)iseed2=iseed2+im2
  j=1+iy/ndiv
  iy=iv(j)-iseed2
  iv(j)=iseed
  if(iy<1)iy=iy+imm1
  ranx=min(am*iy,rnmx)
end function
end module

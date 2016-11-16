 PROGRAM compare_kpoints
!-----------------!
! COMPARE_KPOINTS !
!-----------------! 
 USE constants
 USE utils
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-10
 INTEGER :: i,ierr,size,label,list
 REAL(dp) :: kpoint(3),gvec_frac(3)

 open(unit=7,file='kpoints.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('','Unable to open kpoints.dat file.')
 open(unit=8,file='gvectors_frac.dat',status='old',iostat=ierr)
 if(ierr/=0)call errstop('','Unable to open gvectors_frac.dat file.')
 open(unit=9,file='list.dat',status='new',iostat=ierr)
 if(ierr/=0)call errstop('','Unable to open list.dat file.')

 read(8,*,iostat=ierr)size
 if(ierr/=0)stop

 do i=1,size
  read(8,*,iostat=ierr)label,gvec_frac(1:3)
  if(ierr/=0)call errstop('','Problem reading gvectors_frac.dat file.')
  gvec_frac(1:3)=modulo(0.5d0+gvec_frac(1:3)+tol,1.d0)-0.5d0-tol
  do
   read(7,*,iostat=ierr)list,kpoint(1:3)
   if(ierr>0)then
    call errstop('','Problem reading kpoints.dat file.')
   elseif(ierr<0)then
    exit
   endif ! ierr
   kpoint(1:3)=modulo(0.5d0+kpoint(1:3)+tol,1.d0)-0.5d0-tol
   if(all(abs(gvec_frac(1:3)-kpoint(1:3))<tol))write(9,*)list,label
  enddo
  rewind(7)
 enddo

 close(7) ; close(8) ; close(9)

 END PROGRAM compare_kpoints

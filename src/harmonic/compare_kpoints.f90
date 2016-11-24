module compare_kpoints_module
  implicit none
contains

subroutine compare_kpoints()
  use constants
  use utils
  implicit none
  
  real(dp),parameter :: tol=1.d-10
  integer :: i,ierr,size,label,list
  real(dp) :: kpoint(3),gvec_frac(3)

  open(unit=7,file='kpoints.dat',status='old',iostat=ierr)
  if (ierr/=0) call errstop('','unable to open kpoints.dat file.')
  open(unit=8,file='gvectors_frac.dat',status='old',iostat=ierr)
  if (ierr/=0) call errstop('','unable to open gvectors_frac.dat file.')
  open(unit=9,file='list.dat',status='new',iostat=ierr)
  if (ierr/=0) call errstop('','unable to open list.dat file.')

  read(8,*,iostat=ierr) size
  if (ierr/=0) stop

  do i=1,size
    read(8,*,iostat=ierr) label, gvec_frac(1:3)
    if(ierr/=0) call errstop('','problem reading gvectors_frac.dat file.')
    gvec_frac(1:3)=modulo(0.5d0+gvec_frac(1:3)+tol,1.d0)-0.5d0-tol
    do while (ierr==0)
      read(7,*,iostat=ierr)list,kpoint(1:3)
      if(ierr>0)then
        call errstop('','problem reading kpoints.dat file.')
      elseif (ierr==0)then
        kpoint(1:3)=modulo(0.5d0+kpoint(1:3)+tol,1.d0)-0.5d0-tol
        if (all(abs(gvec_frac(1:3)-kpoint(1:3))<tol)) then
          write(9,*) list, label
        endif
      endif ! ierr
    enddo
    rewind(7)
  enddo

  close(7) ; close(8) ; close(9)

  end subroutine
end module

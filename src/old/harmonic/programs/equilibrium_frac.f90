 MODULE constants
!--------------------------------------------------------------!
! Numerical constants and constants for variable declarations. !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 END MODULE constants

 MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 USE constants
 IMPLICIT NONE

 CONTAINS

  SUBROUTINE inv_33(A,B)
    ! This subroutine calculates the inverse B of matrix A.
    ! A and B are real, 3x3 matrices.
    IMPLICIT NONE
    REAL(dp),INTENT(in) :: A(3,3)
    REAL(dp),INTENT(out) :: B(3,3)
    REAL(dp) :: d
    d=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+ &
      &A(2,1)*(A(3,2)*A(1,3)-A(1,2)*A(3,3))+ &
      &A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    IF(d==0.d0)THEN
      WRITE(*,*)'Error in inv_33: singular matrix.'
      STOP
    ENDIF
    d=1.d0/d
    B(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))*d
    B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))*d
    B(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))*d
    B(2,1)=(A(3,1)*A(2,3)-A(2,1)*A(3,3))*d
    B(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))*d
    B(2,3)=(A(2,1)*A(1,3)-A(1,1)*A(2,3))*d
    B(3,1)=(A(2,1)*A(3,2)-A(2,2)*A(3,1))*d
    B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))*d
    B(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))*d
  END SUBROUTINE inv_33

END MODULE utils


program equilibrium_frac
  use constants
  use utils
  implicit none
  integer :: i,no_atoms
  real(dp) :: lattice(3,3),trans_lattice(3,3),inv_lattice(3,3)
  real(dp),allocatable :: atoms(:,:),mass(:),frac_atoms(:,:)
  character(2),allocatable :: species(:)

  ! Read in atomic positions
  open(1,file='super_equilibrium.dat')
  read(1,*)no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms),species(no_atoms),frac_atoms(no_atoms,3))
  do i=1,no_atoms
    read(1,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  close(1)

  ! Read in superlattice
  open(1,file='super_lattice.dat')
  do i=1,3
    read(1,*)lattice(i,:)
  enddo ! i
  close(1) 
  trans_lattice=transpose(lattice)
  call inv_33(trans_lattice,inv_lattice)
  do i=1,no_atoms
    frac_atoms(i,1:3)=atoms(i,1)*inv_lattice(1:3,1)+atoms(i,2)*inv_lattice(1:3,2)+&
                &atoms(i,3)*inv_lattice(1:3,3)
  enddo ! i

  ! Write out fractional atomic positions
  open(1,file='super_equilibrium_frac.dat')
  write(1,*)no_atoms
  do i=1,no_atoms
    write(1,*)species(i),mass(i),frac_atoms(i,:)
  enddo ! i
  close(1)

end program equilibrium_frac

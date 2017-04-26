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

 FUNCTION determinant33(A)
!-----------------------------------------------------!
! Given a 3x3 matrix A, this function returns det(A). !
!-----------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: A(3,3)
 INTEGER :: determinant33

 determinant33=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
  &+A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
  &+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

 END FUNCTION determinant33

 FUNCTION determinant33_real(A)
!-----------------------------------------------------!
! Given a 3x3 matrix A, this function returns det(A). !
!-----------------------------------------------------!
 IMPLICIT NONE
 real(dp),INTENT(in) :: A(3,3)
 real(dp) :: determinant33_real

 determinant33_real=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
  &+A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
  &+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

 END FUNCTION determinant33_real



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

program generate_supercell_kpoint_mesh_qe
  use constants
  use utils
  implicit none
  ! Working variables
  integer :: i,j,k
  real,parameter :: pi=3.14159265358979324d0
  ! Input variables
  integer :: prim_mesh(3),sc_mesh(3),supercell(3,3)
  real(dp),allocatable :: directed_mesh1(:,:),directed_mesh2(:,:),directed_mesh3(:,:)
  real(dp),allocatable :: mesh(:,:)
  real(dp),allocatable :: sc_directed_mesh1(:,:),sc_directed_mesh2(:,:),sc_directed_mesh3(:,:)
  real(dp) :: rec_distance1(3),rec_distance2(3),rec_distance3(3)
  real(dp) :: sc_rec_distance1(3),sc_rec_distance2(3),sc_rec_distance3(3)
  real(dp) :: sc_dist(3),super_lattice(3,3),rec_super_lattice(3,3)
  real(dp) :: dist(3),lattice(3,3),rec_lattice(3,3)

  ! Read in mesh of primitive cell
  open(1,file='kpoints.in')
  read(1,*)
  read(1,*)prim_mesh(:)
  close(1)

  ! Read in primitive lattice
  open(1,file='lattice.dat')
  do i=1,3
    read(1,*)lattice(i,:)
  enddo ! i
  close(1)

  ! Construct reciprocal primitive lattice
  call inv_33(lattice,rec_lattice)
  rec_lattice=2.d0*pi*transpose(rec_lattice)
  do i=1,3
    dist(i)=sqrt(dot_product(rec_lattice(i,:),rec_lattice(i,:)))
  enddo ! i

  ! Read in SC lattice
  open(1,file='super_lattice.dat')
  do i=1,3
    read(1,*)super_lattice(i,:)
  enddo ! i
  close(1)

  ! Construct reciprocal SC lattice
  call inv_33(super_lattice,rec_super_lattice)
  rec_super_lattice=2.d0*pi*transpose(rec_super_lattice)
  do i=1,3
    sc_dist(i)=sqrt(dot_product(rec_super_lattice(i,:),rec_super_lattice(i,:)))
  enddo ! i

  open(1,file='sc_kpoints.dat')
  write(1,*)int(prim_mesh(1)*sc_dist(1)/dist(1))+1,int(prim_mesh(2)*sc_dist(2)/dist(2))+1,&
          &int(prim_mesh(3)*sc_dist(3)/dist(3))+1,0,0,0
  close(1)

end program generate_supercell_kpoint_mesh_qe

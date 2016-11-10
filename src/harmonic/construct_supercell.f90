! Program to construct supercell from primitive cell

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

 LOGICAL FUNCTION is_in_supercell(pos,super_lattice)
   ! This function returns T if and only if pos is in the supercell
   IMPLICIT NONE
   REAL(dp),INTENT(in) :: pos(3),super_lattice(3,3)
   REAL(dp) :: trans_super_lattice(3,3),inv_super_lattice(3,3),frac_pos(3)
   REAL(dp) :: t,a,b,c,f1,f2,f3
   REAL(dp),PARAMETER :: tol=1.d-2

   trans_super_lattice=transpose(super_lattice)
   CALL inv_33(trans_super_lattice,inv_super_lattice)
   frac_pos(1:3)=pos(1)*inv_super_lattice(1:3,1)+pos(2)*inv_super_lattice(1:3,2)+&
                &pos(3)*inv_super_lattice(1:3,3)
   !write(*,*)frac_pos
   if(frac_pos(1)>-tol.and.frac_pos(1)<(1.d0-tol))then
     if(frac_pos(2)>-tol.and.frac_pos(2)<(1.d0-tol))then
       is_in_supercell=(frac_pos(3)>-tol.and.frac_pos(3)<(1.d0-tol))
     else
       is_in_supercell=.FALSE.
     endif ! 2nd component inside
   else
     is_in_supercell=.FALSE.
   endif ! 1st component inside

 END FUNCTION is_in_supercell


END MODULE utils



PROGRAM construct_supercell
  USE utils
  USE constants
  IMPLICIT NONE
  ! Working variables
  INTEGER :: i,j,k,atom_counter,delta
  INTEGER :: dira,dirb,dirc
  REAL(dp) :: pos(3)
  ! Input variables
  INTEGER :: no_atoms,supercell(3,3),sc_size
  REAL(dp) :: lattice(3,3),super_lattice(3,3)
  REAL(dp),ALLOCATABLE :: atoms(:,:),mass(:),super_atoms(:,:),super_mass(:)
  CHARACTER(2),ALLOCATABLE :: species(:),super_species(:)
  

  ! Read in lattice
  open(1,file='lattice.dat')
  do i=1,3 
    read(1,*)lattice(i,1:3)
  enddo ! i
  close(1) ! lattice.dat

  ! Read in equilibrium
  open(1,file='equilibrium.dat')
  read(1,*)no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms))
  allocate(species(no_atoms))
  do i=1,no_atoms
    read(1,*)species(i),mass(i),atoms(i,1:3)
  enddo ! i
  close(1) ! equilibrium.dat


  ! Read in supercell matrix
  open(1,file='supercell.dat')
  do i=1,3
    read(1,*)supercell(i,1:3)
  enddo ! i
  close(1) ! supercell.dat
  sc_size=abs(determinant33(supercell))
  allocate(super_atoms(no_atoms*sc_size,3),super_mass(no_atoms*sc_size),&  
    &super_species(no_atoms*sc_size))

  ! Generate supercell lattice
  super_lattice=MATMUL(supercell,lattice)
  open(1,file='super_lattice.dat')
  do i=1,3
    write(1,*)super_lattice(i,1:3)
  enddo  
  close(1) ! super_lattice

  ! Generate supercell atoms
  atom_counter=0
  delta=sc_size*3
  do i=1,no_atoms
    do dira=-delta,delta
      do dirb=-delta,delta
        do dirc=-delta,delta
          pos(:)=atoms(i,:)+dira*lattice(1,:)+dirb*lattice(2,:)+dirc*lattice(3,:)
          if(is_in_supercell(pos,super_lattice))then
            atom_counter=atom_counter+1
            super_atoms(atom_counter,:)=pos(:)
            super_mass(atom_counter)=mass(i)
            super_species(atom_counter)=species(i)
          endif ! is in supercell          
        enddo ! dirc
      enddo ! dirb
    enddo ! dira
  enddo ! Loop over primitive cell atoms

  if(atom_counter/=(no_atoms*sc_size))then
    write(*,*)'Error in placing atoms in supercell! Please, try increasing delta.'
    STOP
  endif ! error in generating supercell 

  open(1,file='super_equilibrium.dat')
  write(1,*)no_atoms*sc_size
  do i=1,no_atoms*sc_size
    write(1,*)super_species(i),super_mass(i),super_atoms(i,:)
  enddo ! i
  close(1)

END PROGRAM construct_supercell

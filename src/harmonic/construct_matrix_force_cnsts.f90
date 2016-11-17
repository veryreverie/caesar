! This program determines the atomic displacements required to construct the matrix of 
! force constants

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module construct_matrix_force_cnsts_module
  implicit none
contains

subroutine construct_matrix_force_cnsts()
  use utils
  use constants
  use linear_algebra, only : inv_33
  implicit none
  ! Working variables
  integer :: i,j,k
  real(dp) :: rot_pos(3),reduced_pos(3),rot_pos_frac(3)
  logical :: related,yrelated,zrelated,false
  real(dp) :: x(3),y(3),z(3),roty(3),rotz(3)
  real(dp) :: tol=1.d-5,tol_mod=1.d-8
  ! Input variables
  integer :: no_symm,no_atoms
  real(dp),allocatable :: rotation(:,:,:),offset(:,:),atom_pos(:,:),atom_pos_frac(:,:),mass(:)
  real(dp) :: lattice(3,3),inv_lattice(3,3),trans_lattice(3,3),supercell(3,3)
  character(2) :: dump_ch

  false=.false.

  ! Prepare output file
  open(71,file='force_constants.dat')

  ! Read in symmetry operations
  open(1,file='symmetry.dat')
  read(1,*)no_symm
  allocate(rotation(no_symm,3,3))
  allocate(offset(no_symm,3))
  do i=1,no_symm
    read(1,*)rotation(i,1,:)
    read(1,*)rotation(i,2,:)
    read(1,*)rotation(i,3,:)
    read(1,*)offset(i,:)
  enddo ! i
  close(1)

  ! Read in lattice
  open(1,file='lattice.dat')
  read(1,*)lattice(1,:)
  read(1,*)lattice(2,:)
  read(1,*)lattice(3,:)
  close(1)
  trans_lattice=transpose(lattice)
  call inv_33(trans_lattice,inv_lattice)

  ! Read in supercell matrix
  open(1,file='supercell.dat')
  read(1,*)supercell(1,:)
  read(1,*)supercell(2,:)
  read(1,*)supercell(3,:)
  close(1)
  supercell=transpose(supercell)
  call inv_33(supercell,supercell)
  ! Transform offsets to primitive cell coordinates
  do i=1,no_symm
    offset(i,:)=matmul(supercell,offset(i,:))  
  enddo ! i

  ! Read in atomic positions
  open(1,file='super_equilibrium.dat')
  read(1,*)no_atoms
  allocate(atom_pos(no_atoms,3),mass(no_atoms))
  allocate(atom_pos_frac(no_atoms,3))
  do i=1,no_atoms
    read(1,*)dump_ch,mass(i),atom_pos(i,:)
  enddo ! i
  close(1)

  do i=1,no_atoms
    atom_pos_frac(i,:)=atom_pos(i,1)*inv_lattice(:,1)+atom_pos(i,2)*inv_lattice(:,2)+&
                    &atom_pos(i,3)*inv_lattice(:,3)
  enddo ! i

  ! Check which Cartesian directions are needed
  yrelated=.false.
  zrelated=.false.
  x(1)=1; x(2)=0; x(3)=0
  y(1)=0; y(2)=1; y(3)=0
  z(1)=0; z(2)=0; z(3)=1
  do i=1,no_symm
    roty=matmul(rotation(i,:,:),y)
    rotz=matmul(rotation(i,:,:),z)
    if(abs(roty(1)-x(1))<tol.and.abs(roty(2)-x(2))<tol.and.abs(roty(3)-x(3))<tol)yrelated=.true.
    if(abs(rotz(1)-x(1))<tol.and.abs(rotz(2)-x(2))<tol.and.abs(rotz(3)-x(3))<tol)zrelated=.true.
  enddo ! i
  

  ! Apply point group to atomic positions
  do i=1,no_atoms
    related=.false.
    !if(i==2)then
    do j=1,no_symm 
      rot_pos=matmul(rotation(j,:,:),atom_pos(i,:))
      rot_pos_frac(:)=rot_pos(1)*inv_lattice(:,1)+rot_pos(2)*inv_lattice(:,2)+&
                     &rot_pos(3)*inv_lattice(:,3)+offset(j,:)
      !write(*,*)j,rot_pos_frac,atom_pos_frac(1,:)
      do k=1,i-1 
        reduced_pos(:)=modulo(abs(rot_pos_frac(:)-atom_pos_frac(k,:)-tol_mod),1.d0) 
     !   write(*,*)i,j,reduced_pos
        if((abs(reduced_pos(1))<tol.or.abs(reduced_pos(1)-1)<tol.or.abs(reduced_pos(1)+1)<tol)&
     &.and.(abs(reduced_pos(2))<tol.or.abs(reduced_pos(2)-1)<tol.or.abs(reduced_pos(2)+1)<tol)&
     &.and.(abs(reduced_pos(3))<tol.or.abs(reduced_pos(3)-1)<tol.or.abs(reduced_pos(3)+1)<tol))&
     &then
          related=.true.
        endif 
      enddo ! k
    enddo ! j
    !endif ! i==2
    if(related.eqv.false)then
      write(71,*)i,1 
      if(yrelated.eqv.false)write(71,*)i,2
      if(zrelated.eqv.false)write(71,*)i,3
    endif
  enddo ! i

end subroutine
end module

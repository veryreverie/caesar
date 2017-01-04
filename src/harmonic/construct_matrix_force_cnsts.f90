! This program determines the atomic displacements required to construct the
! matrix of force constants

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module construct_matrix_force_cnsts_module
  implicit none
contains

subroutine construct_matrix_force_cnsts(filenames)
  use constants,      only : dp
  use linear_algebra, only : inv_33
  use file_io,        only : open_read_file, open_write_file
  implicit none
  
  character(100), intent(in) :: filenames(:)
  
  ! Input variables
  integer :: no_symm,no_atoms
  real(dp),allocatable :: rotation(:,:,:),offset(:,:),atom_pos(:,:),atom_pos_frac(:,:),mass(:)
  real(dp) :: lattice(3,3),inv_lattice(3,3),trans_lattice(3,3),supercell(3,3)
  character(2) :: dump_ch
  
  ! Parameters
  real(dp), parameter :: x(3) = (/ 1.d0, 0.d0, 0.d0 /)
  real(dp), parameter :: y(3) = (/ 0.d0, 1.d0, 0.d0 /)
  real(dp), parameter :: z(3) = (/ 0.d0, 0.d0, 1.d0 /)
  real(dp), parameter :: tol=1.d-5
  real(dp), parameter :: tol_mod=1.d-8
  
  ! Working variables
  integer :: i,j,k
  real(dp) :: rot_pos(3),reduced_pos(3),rot_pos_frac(3)
  logical :: related,yrelated,zrelated
  real(dp) :: roty(3),rotz(3)
  
  ! file units
  integer :: symmetry_file
  integer :: lattice_file
  integer :: supercell_file
  integer :: super_equilibrium_file
  integer :: force_constants_file
  
  ! Read in symmetry operations
  symmetry_file = open_read_file(filenames(1))
  read(symmetry_file,*)no_symm
  allocate(rotation(no_symm,3,3))
  allocate(offset(no_symm,3))
  do i=1,no_symm
    read(symmetry_file,*)rotation(i,1,:)
    read(symmetry_file,*)rotation(i,2,:)
    read(symmetry_file,*)rotation(i,3,:)
    read(symmetry_file,*)offset(i,:)
  enddo ! i
  close(symmetry_file)

  ! Read in lattice
  lattice_file = open_read_file(filenames(2))
  read(lattice_file,*)lattice(1,:)
  read(lattice_file,*)lattice(2,:)
  read(lattice_file,*)lattice(3,:)
  close(lattice_file)
  
  trans_lattice=transpose(lattice)
  call inv_33(trans_lattice,inv_lattice)

  ! Read in supercell matrix
  supercell_file = open_read_file(filenames(3))
  read(supercell_file,*)supercell(1,:)
  read(supercell_file,*)supercell(2,:)
  read(supercell_file,*)supercell(3,:)
  close(supercell_file)
  
  supercell=transpose(supercell)
  call inv_33(supercell,supercell)
  ! Transform offsets to primitive cell coordinates
  do i=1,no_symm
    offset(i,:)=matmul(supercell,offset(i,:))  
  enddo ! i

  ! Read in atomic positions
  super_equilibrium_file = open_read_file(filenames(4))
  read(super_equilibrium_file,*)no_atoms
  allocate(atom_pos(no_atoms,3),mass(no_atoms))
  allocate(atom_pos_frac(no_atoms,3))
  do i=1,no_atoms
    read(super_equilibrium_file,*)dump_ch,mass(i),atom_pos(i,:)
  enddo ! i
  close(super_equilibrium_file)

  do i=1,no_atoms
    atom_pos_frac(i,:) = atom_pos(i,1)*inv_lattice(:,1) &
                     & + atom_pos(i,2)*inv_lattice(:,2) &
                     & + atom_pos(i,3)*inv_lattice(:,3)
  enddo ! i

  ! Check which Cartesian directions are needed
  yrelated=.false.
  zrelated=.false.
  do i=1,no_symm
    roty=matmul(rotation(i,:,:),y)
    rotz=matmul(rotation(i,:,:),z)
    if( abs(roty(1)-x(1))<tol .and. &
      & abs(roty(2)-x(2))<tol .and. &
      & abs(roty(3)-x(3))<tol ) yrelated=.true.
    if( abs(rotz(1)-x(1))<tol .and. &
      & abs(rotz(2)-x(2))<tol .and. &
      & abs(rotz(3)-x(3))<tol ) zrelated=.true.
  enddo ! i
  
  ! Prepare output file
  force_constants_file = open_write_file(filenames(5))

  ! Apply point group to atomic positions
  do i=1,no_atoms
    related=.false.
    do j=1,no_symm 
      rot_pos=matmul(rotation(j,:,:),atom_pos(i,:))
      rot_pos_frac(:) = rot_pos(1)*inv_lattice(:,1) &
                    & + rot_pos(2)*inv_lattice(:,2) &
                    & + rot_pos(3)*inv_lattice(:,3) &
                    & + offset(j,:)
      do k=1,i-1
        reduced_pos = rot_pos_frac(:)-atom_pos_frac(k,:)-tol_mod
        reduced_pos = dabs(reduced_pos-nint(reduced_pos))
        if ( reduced_pos(1)<tol .and. &
           & reduced_pos(2)<tol .and. &
           & reduced_pos(3)<tol) then
          related=.true.
          exit
        endif 
      enddo ! k
      if (related) exit
    enddo ! j
    
    if (.not. related) then
      write(force_constants_file,*) i, 1 
      if (.not. yrelated) write(force_constants_file,*) i, 2
      if (.not. zrelated) write(force_constants_file,*) i, 3
    endif
  enddo ! i
  
  close(force_constants_file)

end subroutine
end module

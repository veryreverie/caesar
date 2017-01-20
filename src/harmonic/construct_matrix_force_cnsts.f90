! This program determines the atomic displacements required to construct the
! matrix of force constants

! The symmetry operations to consider are nicely outlined here:
! http://www.homepages.ucl.ac.uk/~ucfbdxa/phon/node4.html

module construct_matrix_force_cnsts_module
  implicit none
contains

subroutine construct_matrix_force_cnsts(filenames)
  use constants,        only : dp
  use linear_algebra,   only : inv_33
  use file_io,          only : open_read_file, open_write_file
  use structure_module, only : StructureData, read_structure_file, drop
  use string_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  ! Input variables
  integer :: no_symm
  real(dp),allocatable :: rotation(:,:,:),offset(:,:),atom_pos_frac(:,:)
  real(dp) :: supercell(3,3)
  
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
  
  type(StructureData) :: structure
  type(StructureData) :: superstructure
  
  ! file units
  integer :: symmetry_file
  integer :: supercell_file
  integer :: force_constants_file
  
  ! filenames
  character(100) :: symmetry_filename
  character(100) :: structure_filename
  character(100) :: supercell_filename
  character(100) :: superstructure_filename
  character(100) :: force_constants_filename
  
  ! Read filenames from input
  symmetry_filename = filenames(1)
  structure_filename = filenames(2)
  supercell_filename = filenames(3)
  superstructure_filename = filenames(4)
  force_constants_filename = filenames(5)
  
  ! Read in symmetry operations
  symmetry_file = open_read_file(symmetry_filename)
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
  structure = read_structure_file(structure_filename)

  ! Read in supercell matrix
  supercell_file = open_read_file(supercell_filename)
  read(supercell_file,*)supercell(1,:)
  read(supercell_file,*)supercell(2,:)
  read(supercell_file,*)supercell(3,:)
  close(supercell_file)
  
  supercell = inv_33(transpose(supercell))
  ! Transform offsets to primitive cell coordinates
  do i=1,no_symm
    offset(i,:)=matmul(supercell,offset(i,:))  
  enddo ! i

  ! Read in atomic positions
  superstructure = read_structure_file(superstructure_filename)

  atom_pos_frac = matmul(structure%recip_lattice,superstructure%atoms)

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
  force_constants_file = open_write_file(force_constants_filename)

  ! Apply point group to atomic positions
  do i=1,superstructure%no_atoms
    related=.false.
    do j=1,no_symm 
      rot_pos=matmul(rotation(j,:,:),superstructure%atoms(:,i))
      rot_pos_frac(:) = rot_pos(1)*structure%recip_lattice(:,1) &
                    & + rot_pos(2)*structure%recip_lattice(:,2) &
                    & + rot_pos(3)*structure%recip_lattice(:,3) &
                    & + offset(j,:)
      do k=1,i-1
        reduced_pos = rot_pos_frac(:)-atom_pos_frac(:,k)-tol_mod
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
  
  call drop(structure)
end subroutine
end module

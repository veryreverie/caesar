module generate_supercell_kpoint_mesh_qe_module
  implicit none
contains

subroutine generate_supercell_kpoint_mesh_qe(filenames)
  use constants
  use utils
  use linear_algebra, only : inv_33
  use file_io,        only : open_read_file, open_write_file
  implicit none
  
  character(100), intent(in) :: filenames(:)
  
  ! Working variables
  integer :: i,j,k
!  real,parameter :: pi=3.14159265358979324d0
  ! Input variables
  integer :: prim_mesh(3),sc_mesh(3),supercell(3,3)
  real(dp),allocatable :: directed_mesh1(:,:),directed_mesh2(:,:),directed_mesh3(:,:)
  real(dp),allocatable :: mesh(:,:)
  real(dp),allocatable :: sc_directed_mesh1(:,:),sc_directed_mesh2(:,:),sc_directed_mesh3(:,:)
  real(dp) :: rec_distance1(3),rec_distance2(3),rec_distance3(3)
  real(dp) :: sc_rec_distance1(3),sc_rec_distance2(3),sc_rec_distance3(3)
  real(dp) :: sc_dist(3),super_lattice(3,3),rec_super_lattice(3,3)
  real(dp) :: dist(3),lattice(3,3),rec_lattice(3,3)
  
  ! file units
  integer :: kpoints_file
  integer :: lattice_file
  integer :: super_lattice_file
  integer :: sc_kpoints_file
  
  ! Read in mesh of primitive cell
  kpoints_file = open_read_file(filenames(1))
  read(kpoints_file,*)
  read(kpoints_file,*) prim_mesh(:)
  close(kpoints_file)

  ! Read in primitive lattice
  lattice_file = open_read_file(filenames(2))
  do i=1,3
    read(lattice_file,*) lattice(i,:)
  enddo
  close(lattice_file)

  ! Construct reciprocal primitive lattice
  call inv_33(lattice,rec_lattice)
  rec_lattice=2.d0*pi*transpose(rec_lattice)
  do i=1,3
    dist(i)=sqrt(dot_product(rec_lattice(i,:),rec_lattice(i,:)))
  enddo ! i

  ! Read in SC lattice
  super_lattice_file = open_read_file(filenames(3))
  do i=1,3
    read(super_lattice_file,*) super_lattice(i,:)
  enddo
  close(super_lattice_file)

  ! Construct reciprocal SC lattice
  call inv_33(super_lattice,rec_super_lattice)
  rec_super_lattice=2.d0*pi*transpose(rec_super_lattice)
  do i=1,3
    sc_dist(i)=sqrt(dot_product(rec_super_lattice(i,:),rec_super_lattice(i,:)))
  enddo ! i
  
  sc_kpoints_file = open_write_file(filenames(4))
  write(sc_kpoints_file,*) int(prim_mesh(1)*sc_dist(1)/dist(1))+1, &
                         & int(prim_mesh(2)*sc_dist(2)/dist(2))+1, &
                         & int(prim_mesh(3)*sc_dist(3)/dist(3))+1, &
                         & 0,                                      &
                         & 0,                                      &
                         & 0
  close(sc_kpoints_file)
end subroutine
end module

module generate_supercell_kpoint_mesh_qe_module
  implicit none
contains

subroutine generate_supercell_kpoint_mesh_qe(filenames)
  use constants
  use utils
  use linear_algebra,   only : inv_33
  use file_io,          only : open_read_file, open_write_file
  use structure_module, only : StructureData, read_structure_file, drop
  use string_module
  implicit none
  
  type(String), intent(in) :: filenames(:)
  
  ! Working variables
  integer :: i
  
  ! Input variables
  integer :: prim_mesh(3)
  real(dp) :: sc_dist(3),super_lattice(3,3),rec_super_lattice(3,3)
  real(dp) :: dist(3),rec_lattice(3,3)
  
  type(StructureData) :: structure
  
  ! file units
  integer :: kpoints_file
  integer :: super_lattice_file
  integer :: sc_kpoints_file
  
  ! filenames
  character(100) :: kpoints_filename
  character(100) :: structure_filename
  character(100) :: super_lattice_filename
  character(100) :: sc_kpoints_filename
  
  ! The first line of kpoints_filename
  character(100) :: header
  
  ! Read filenames from input
  kpoints_filename = filenames(1)
  structure_filename = filenames(2)
  super_lattice_filename = filenames(3)
  sc_kpoints_filename = filenames(4)
  
  ! Read in mesh of primitive cell
  kpoints_file = open_read_file(kpoints_filename)
  read(kpoints_file,"(a)") header
  read(kpoints_file,*) prim_mesh(:)
  close(kpoints_file)

  ! Construct reciprocal primitive lattice
  rec_lattice = 2.d0*pi*structure%recip_lattice
  do i=1,3
    dist(i) = norm2(rec_lattice(i,:))
  enddo

  ! Read in SC lattice
  super_lattice_file = open_read_file(super_lattice_filename)
  do i=1,3
    read(super_lattice_file,*) super_lattice(i,:)
  enddo
  close(super_lattice_file)

  ! Construct reciprocal SC lattice
  rec_super_lattice=2.d0*pi*transpose(inv_33(super_lattice))
  do i=1,3
    sc_dist(i)=sqrt(dot_product(rec_super_lattice(i,:),rec_super_lattice(i,:)))
  enddo ! i
  
  sc_kpoints_file = open_write_file(sc_kpoints_filename)
  write(sc_kpoints_file,"(a)") trim(header)
  write(sc_kpoints_file,*) int(prim_mesh(1)*sc_dist(1)/dist(1))+1, &
                         & int(prim_mesh(2)*sc_dist(2)/dist(2))+1, &
                         & int(prim_mesh(3)*sc_dist(3)/dist(3))+1, &
                         & 0,                                      &
                         & 0,                                      &
                         & 0
  close(sc_kpoints_file)
  
  call drop(structure)
end subroutine
end module

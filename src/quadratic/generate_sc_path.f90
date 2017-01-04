module generate_sc_path_module
  implicit none
contains

subroutine generate_sc_path(filenames)
  use constants, only : dp
  use file_io,   only : open_read_file, open_write_file
  implicit none
  
  character(32), intent(in) :: filenames(:)
  
  integer :: i,no_points
  integer :: supercell(3,3)
  real(dp),allocatable :: path(:,:),sc_path(:,:)
  
  ! file units
  integer :: supercell_file
  integer :: bs_path_file
  integer :: sc_bs_path_file

  ! Read in supercell matrix
  supercell_file = open_read_file(filenames(1))
  do i=1,3
    read(supercell_file,*)supercell(i,:)
  enddo
  close(supercell_file)

  ! Read in path
  bs_path_file = open_read_file(filenames(2))
  sc_bs_path_file = open_write_file(filenames(3))
  read(bs_path_file,*) no_points
  allocate(path(no_points,3),sc_path(no_points,3))
  do i=1,no_points
    read(bs_path_file,*)path(i,:)
    sc_path(i,:)=matmul(supercell,path(i,:))
    write(sc_bs_path_file,*)sc_path(i,:)
  enddo
  close(bs_path_file)
  close(sc_bs_path_file)
end subroutine
end module

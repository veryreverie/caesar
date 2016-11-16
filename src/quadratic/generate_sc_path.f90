module generate_sc_path_module
  implicit none
contains

subroutine generate_sc_path
  implicit none
  integer :: i,no_points
  integer :: supercell(3,3)
  real,allocatable :: path(:,:),sc_path(:,:)

  ! Read in supercell matrix
  open(1,file='supercell.dat')
  do i=1,3
    read(1,*)supercell(i,:)
  enddo ! i
  close(1)

  ! Read in path
  open(1,file='bs_path.dat')
  open(2,file='sc_bs_path.dat')
  read(1,*)no_points
  allocate(path(no_points,3),sc_path(no_points,3))
  do i=1,no_points
    read(1,*)path(i,:)
    sc_path(i,:)=matmul(supercell,path(i,:))
    write(2,*)sc_path(i,:)
  enddo ! i
  close(1)
  close(2)
end subroutine
end module

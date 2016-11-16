module construct_finite_displacement_module
  implicit none
contains

subroutine construct_finite_displacement()
  implicit none
  integer,parameter :: dp=kind(1.d0)
  ! Working variables
  integer :: i
  ! Input variables
  integer :: atom,disp,no_atoms
  real(dp) :: lattice(3,3)
  real(dp),allocatable :: atoms(:,:),mass(:)
  character(2),allocatable :: species(:)

  ! Read in displacement
  open(1,file='disp.dat')
  read(1,*)atom,disp 
  close(1)

  ! Read in structure
  open(1,file='super_lattice.dat')
  read(1,*)lattice(1,:)
  read(1,*)lattice(2,:)
  read(1,*)lattice(3,:)
  close(1)
  open(1,file='super_equilibrium.dat')
  read(1,*)no_atoms
  allocate(atoms(no_atoms,3),mass(no_atoms),species(no_atoms))
  do i=1,no_atoms 
    read(1,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  close(1)

  ! Generate distorted structure
  open(1,file='positive.dat')
  open(2,file='negative.dat')
  write(1,*)'Lattice'
  write(2,*)'Lattice'
  do i=1,3
    write(1,*)lattice(i,:)
    write(2,*)lattice(i,:)
  enddo ! i
  write(1,*)'Atoms'
  write(2,*)'Atoms'
  do i=1,no_atoms
    if(i==atom)then
      atoms(i,disp)=atoms(i,disp)+0.01d0
    endif ! i==atom
    write(1,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  do i=1,no_atoms
    if(i==atom)then
      atoms(i,disp)=atoms(i,disp)-0.02d0
    endif ! i==atom
    write(2,*)species(i),mass(i),atoms(i,:)
  enddo ! i
  write(1,*)'Symmetry'
  write(2,*)'Symmetry'
  write(1,*)'End'
  write(2,*)'End'
  close(1)
  close(2)

end subroutine
end module

program combine_forces
  implicit none
  integer,parameter :: dp=kind(1.d0)
  real(dp),allocatable :: positive(:),negative(:)
  REAL(dp),PARAMETER :: convert_eV_per_A_to_au=0.01944689814638725057d0
  integer,allocatable :: atom1(:),disp1(:),atom2(:),disp2(:)
  integer :: i,no_atoms
  

  ! Read in number of atoms
  open(1,file='super_equilibrium.dat')
  read(1,*)no_atoms
  close(1)

  ! Read in force constants
  allocate(positive(no_atoms*3),negative(no_atoms*3))
  allocate(atom1(no_atoms*3),disp1(no_atoms*3))
  allocate(atom2(no_atoms*3),disp2(no_atoms*3))
  open(1,file='positive.dat')
  open(2,file='negative.dat')
  do i=1,no_atoms*3
    read(1,*)atom1(i),disp1(i),atom2(i),disp2(i),positive(i)
    read(2,*)atom1(i),disp1(i),atom2(i),disp2(i),negative(i)
  enddo ! i
  close(1)
  close(2)

  ! Write out combined forces constants
  open(1,file='forces.dat')
  do i=1,no_atoms*3
    write(1,*)atom1(i),disp1(i),atom2(i),disp2(i),(positive(i)-negative(i))/(2.d0*0.01d0)*convert_eV_per_A_to_au
  enddo ! i
  close(1)

end program combine_forces

program convert_forces_from_Rybohr_to_eVang
  implicit none
  integer :: i
  integer :: no_atoms,in1,in2,in3,in4
  real :: force


  open(1,file='no_atoms.txt')
  read(1,*)no_atoms
  close(1)

  open(1,file='forces.dat')
  open(2,file='forces_temp.dat')
  do i=1,no_atoms*3
    read(1,*)in1,in2,in3,in4,force
    write(2,*)in1,in2,in3,in4,force*13.605698066/0.529177
  enddo ! i
  close(1)
  close(2)

end program

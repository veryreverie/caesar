program calculate_gap
  implicit none
  integer :: i,no_points
  character(80) :: filename1,filename2
  real,allocatable :: temp(:),top(:),bottom(:)
 
  write(*,*)'What is the file of the top band?'
  read(*,*)filename1
  write(*,*)'What is the file of the bottom band?'
  read(*,*)filename2
  write(*,*)'How many data points are there?'
  read(*,*)no_points
  allocate(temp(no_points))
  allocate(top(no_points))
  allocate(bottom(no_points))

  open(1,file=filename1)
  open(2,file=filename2)
  open(3,file='gap.dat')
  do i=1,no_points
    read(1,*)temp(i),top(i)
    read(2,*)temp(i),bottom(i)
    write(3,*)temp(i),top(i)-bottom(i)
  enddo ! i
  close(1)
  close(2)
  close(3)
  
  
end program calculate_gap

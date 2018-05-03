real function quad_amplitude(frequency)
    ! The normal mode amplitudes to sample are calculated 
    implicit none
    real,parameter :: tolerance=1.d-6
    real :: frequency, temperature, thermal_energy

     quad_amplitude=sqrt(0.5/frequency)

end function quad_amplitude

! Program to generate amplitudes for energy files
program generate_amplitudes
  implicit none
  ! Input variables
  real :: max_amplitude,frequency
  real,allocatable :: energy(:)
  integer :: first_amplitude,last_amplitude,no_amplitudes,cell_size,mid_amplitude
  ! Working variables
  integer :: i
  real :: amplitude,damplitude,quad_amplitude

  open(1,file='mapping.dat')
  read(1,*)max_amplitude
  read(1,*)first_amplitude,last_amplitude
  close(1)
  no_amplitudes=last_amplitude-first_amplitude+1
  mid_amplitude=(no_amplitudes-1)/2+1
  allocate(energy(no_amplitudes))

  open(1,file='working_frequency.dat')
  read(1,*)frequency
  close(1)
  frequency=frequency/27.21138602

  open(1,file='working_energy.dat')
  do i=1,no_amplitudes
    read(1,*)energy(i)
  enddo ! i
  close(1)
 
  open(1,file='working_size.dat')
  read(1,*)cell_size
  close(1)

  open(1,file='amplitude_energy.dat')
  amplitude=-max_amplitude*quad_amplitude(abs(frequency))
  damplitude=abs(amplitude/first_amplitude)
  do i=1,no_amplitudes
    write(1,*)amplitude,(energy(i)-energy(mid_amplitude))/cell_size
    amplitude=amplitude+damplitude
  enddo ! i
  close(1)


end program generate_amplitudes

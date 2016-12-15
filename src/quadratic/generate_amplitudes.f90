module generate_amplitudes_module
  implicit none
contains

! Program to generate amplitudes for energy files
subroutine generate_amplitudes()
  use constants, only : dp
  implicit none
  
  ! Input variables
  real(dp) :: max_amplitude,frequency
  real(dp),allocatable :: energy(:)
  integer :: first_amplitude,last_amplitude,no_amplitudes,cell_size,mid_amplitude
  ! Working variables
  integer :: i
  real(dp) :: amplitude,damplitude

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
  amplitude=-max_amplitude*sqrt(0.5d0/dabs(frequency))
  damplitude=dabs(amplitude/first_amplitude)
  do i=1,no_amplitudes
    write(1,*)amplitude,(energy(i)-energy(mid_amplitude))/cell_size
    amplitude=amplitude+damplitude
  enddo ! i
  close(1)
end subroutine

! as above, but takes arguments instead of reading from files
pure function generate_amplitudes2(map,energy,frequency,cell_size)&
    & result(output)
  use constants,      only : dp, eV
  use mapping_module, only : Mapping
  implicit none
  
  type(Mapping), intent(in) :: map
  real(dp),      intent(in) :: energy(:)
  real(dp),      intent(in) :: frequency
  integer,       intent(in) :: cell_size
  real(dp), allocatable     :: output(:,:)
  
  ! Working variables
  integer :: i
  real(dp) :: amplitude
  real(dp) :: damplitude
  
  allocate(output(2,map%count))
  amplitude=-map%max*sqrt(0.5d0/dabs(frequency/eV))
  damplitude=dabs(amplitude/map%first)
  do i=1,map%count
    output(1,i) = amplitude
    output(2,i) = (energy(i)-energy(map%mid))/cell_size
    amplitude=amplitude+damplitude
  enddo ! i
end function

end module

module generate_amplitudes_module
  implicit none
contains

! Program to generate amplitudes for energy files
pure function generate_amplitudes(mapping,energy,frequency,cell_size)&
  & result(output)
  use constants, only : dp, eV
  use mapping_module
  implicit none
  
  type(MappingData), intent(in) :: mapping
  real(dp),          intent(in) :: energy(:)
  real(dp),          intent(in) :: frequency
  integer,           intent(in) :: cell_size
  real(dp), allocatable         :: output(:,:)
  
  ! Working variables
  integer  :: i
  real(dp) :: amplitude
  
  allocate(output(2,mapping%count))
  amplitude = -(mapping%max/2.0d0)*(ev/dabs(frequency))
  do i=1,mapping%count
    output(1,i) = amplitude+(i-1)*dabs(amplitude/mapping%first)
    output(2,i) = (energy(i)-energy(mapping%mid))/cell_size
  enddo ! i
end function

end module
